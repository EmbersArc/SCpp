
#include "active_model.hpp"
#include "EcosWrapper.hpp"

#include <iostream>
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using std::array;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::setw;
using std::setfill;

string get_output_path() {
    return "../output/" + Model::get_name() + "/";
}

void make_output_path() {
    string path = "mkdir -p " + get_output_path();
    system(path.c_str());
}


int main() {
    make_output_path();
    Model model;

    const double dt = 1 / double(K-1);

    const double weight_trust_region_sigma = 1;
    const double weight_trust_region_xu = 1;
    const double weight_virtual_control = 100;

    const size_t n_states = Model::n_states;
    const size_t n_inputs = Model::n_inputs;

    Eigen::Matrix<double, n_states, K> X;
    Eigen::Matrix<double, n_inputs, K> U;
    model.initialize(X, U);
    X.setRandom();
    U.setRandom();

    double sigma = model.total_time_guess();


    array<Model::StateMatrix,   (K-1)> A_bar;
    array<Model::ControlMatrix, (K-1)> B_bar;
    array<Model::ControlMatrix, (K-1)> C_bar;
    array<Model::StateVector,   (K-1)> Sigma_bar;
    array<Model::StateVector,   (K-1)> z_bar;



    /** Solver setup **/
    EcosWrapper solver;
    {
        solver.create_tensor_variable("X", {n_states, K}); // states
        solver.create_tensor_variable("U", {n_inputs, K}); // inputs
        solver.create_tensor_variable("sigma", {}); // total time
        solver.create_tensor_variable("Delta_sigma", {}); // squared change of sigma
        solver.create_tensor_variable("Delta", {K}); // squared change of the stacked [ x(k), u(k) ] vector
        solver.create_tensor_variable("norm2_Delta", {}); // 2-norm of the Delta(k) variables

        // shortcuts to access solver variables and create parameters
        auto var = [&](const string &name, const vector<size_t> &indices){ return solver.get_variable(name,indices); };
        auto param = [](double &param_value){ return optimization_problem::Parameter(&param_value); };
        auto param_fn = [](std::function<double()> callback){ return optimization_problem::Parameter(callback); };




        // Build sigma trust region
        // (sigma-sigma0) * (sigma-sigma0) <= Delta_sigma
        //          is equivalent to
        //   norm2([
        //       ((-sigma0)*sigma   +(-0.5)*Delta_sigma  +(0.5+0.5*sigma0*sigma0)  ),
        //       sigma
        //   ])
        //   <= (  +(sigma0)*sigma  +(0.5)*Delta_sigma    +(0.5-0.5*sigma0*sigma0)   )
        {
            // Formulas involving sigma, from the above comment
            auto sigma_fn1 = param_fn([&sigma](){ return -sigma; });
            auto sigma_fn2 = param_fn([&sigma](){ return (0.5+0.5*sigma*sigma); });
            auto sigma_fn3 = param_fn([&sigma](){ return sigma; });
            auto sigma_fn4 = param_fn([&sigma](){ return (0.5-0.5*sigma*sigma); });

            solver.add_constraint( 
                optimization_problem::norm2({
                    sigma_fn1*var("sigma", {})  +  (-0.5)*var("Delta_sigma", {})  +  sigma_fn2,
                    (1.0)*var("sigma", {})
                }) 
                <= sigma_fn3*var("sigma", {})  +  (0.5)*var("Delta_sigma", {})  +  sigma_fn4
            );

            // Minimize Delta_sigma
            solver.add_minimization_term( weight_trust_region_sigma * var("Delta_sigma", {}) );
        }



        
        for (size_t k = 0; k < K; k++) {
            /* 
             * Build state and input trust-region:
             *     (x - x0)^T * (x - x0)  +  (u - u0)^T * (u - u0)  <=  Delta
             * the index k is omitted, but applies to all terms in the constraint.
             * The constraint is equivalent to the SOCP form:
             * 
             * norm2(
             *         ( (-x0^T)*x  +(-u0^T)*u  +(-0.5)*Delta  +(0.5 + 0.5*x0^T*x0 + 0.5*u0^T*u0)),
             *         (I)*x, 
             *         (I)*u
             * )
             *     <= (  ( x0^T)*x  +( u0^T)*u  +( 0.5)*Delta  +(0.5 - 0.5*x0^T*x0 - 0.5*u0^T*u0));
             * 
             */

            vector<optimization_problem::AffineExpression> norm2_args;

            { // (-x0^T)*x  +(-u0^T)*u  +(-0.5)*Delta  +(0.5 + 0.5*x0^T*x0 + 0.5*u0^T*u0)
                optimization_problem::AffineExpression norm2_first_arg;

                // (-x0^T)*x
                for (size_t i = 0; i < n_states; ++i) {
                    norm2_first_arg = norm2_first_arg + 
                        param_fn([&X,i,k](){ return -X(i,k); }) * var("X", {i,k});
                }

                // +(-u0^T)*u
                for (size_t i = 0; i < n_inputs; ++i) {
                    norm2_first_arg = norm2_first_arg + 
                        param_fn([&U,i,k](){ return -U(i,k); }) * var("U", {i,k});
                }

                // +(-0.5)*Delta
                norm2_first_arg = norm2_first_arg + (-0.5) * var("Delta", {k});

                // +(0.5 + 0.5*x0^T*x0 + 0.5*u0^T*u0)
                norm2_first_arg = norm2_first_arg + param_fn([&X,&U,k](){
                    return 0.5 * (1.0 + X.col(k).dot(X.col(k)) + U.col(k).dot(U.col(k)));
                });

                norm2_args.push_back(norm2_first_arg);
            }

            // (I)*x,
            for (size_t i = 0; i < n_states; ++i) {
                norm2_args.push_back( (1.0) * var("X", {i,k}) );
            }

            // (I)*u
            for (size_t i = 0; i < n_inputs; ++i) {
                norm2_args.push_back( (1.0) * var("U", {i,k}) );
            }

            // Right hand side
            // ( x0^T)*x  +( u0^T)*u  +( 0.5)*Delta  +(0.5 - 0.5*x0^T*x0 - 0.5*u0^T*u0)
            optimization_problem::AffineExpression rhs;

            // ( x0^T)*x
            for (size_t i = 0; i < n_states; ++i) {
                rhs = rhs + param(X(i,k)) * var("X", {i,k});
            }

            // +( u0^T)*u
            for (size_t i = 0; i < n_inputs; ++i) {
                rhs = rhs + param(U(i,k)) * var("U", {i,k});
            }

            // +( 0.5)*Delta
            rhs = rhs + (0.5) * var("Delta", {k});

            // +(0.5 - 0.5*x0^T*x0 - 0.5*u0^T*u0)
            rhs = rhs + param_fn([&X,&U,k](){
                return 0.5 * (1.0 - X.col(k).dot(X.col(k)) - U.col(k).dot(U.col(k)));
            });

            solver.add_constraint( optimization_problem::norm2(norm2_args) <= rhs );
        }

        /*
         * Build combined state/input trust region over all K:
         *   norm2([ Delta(1), Delta(2), ... , Delta(K) ]) <= norm2_Delta
         */
        {
            vector<optimization_problem::AffineExpression> norm2_args;
            for (size_t k = 0; k < K; k++) {
                norm2_args.push_back( (1.0) * var("Delta", {k}) );
            }
            solver.add_constraint( optimization_problem::norm2(norm2_args) <= (1.0) * var("norm2_Delta", {}) );

            // Minimize norm2_Delta
            solver.add_minimization_term( weight_trust_region_xu * var("norm2_Delta", {}) );
        }




        model.add_application_constraints(solver);
        solver.compile_problem_structure();
    } /** End solver setup **/


    // Cache indices for performance
    const size_t sigma_index = solver.get_tensor_variable_index("sigma", {});
    size_t X_indices[n_states][K];
    size_t U_indices[n_inputs][K];
    for (size_t k = 0; k < K; k++) {
        for (size_t i = 0; i < n_states; ++i) X_indices[i][k] = solver.get_tensor_variable_index("X",{i,k});
        for (size_t i = 0; i < n_inputs; ++i) U_indices[i][k] = solver.get_tensor_variable_index("U",{i,k});
    }


    const size_t iterations = 40;
    for(size_t it = 0; it < iterations; it++) {
        cout << "Iteration " << it << endl;
        // Write problem to file

        string file_name_prefix;
        {
            ostringstream file_name_prefix_ss;
            file_name_prefix_ss << get_output_path() << "iteration"
            << setfill('0') << setw(3) << it << "_";
            file_name_prefix = file_name_prefix_ss.str();
        }
        
        {
            ofstream f(file_name_prefix + "problem.txt");
            solver.print_problem(f);
        }

        /************************************************************************************/
        cout << "Solving problem" << endl;
        solver.solve_problem();

        // Read solution
        for (size_t k = 0; k < K; k++) {
            for (size_t i = 0; i < n_states; ++i) X(i,k) = solver.get_solution_value(X_indices[i][k]);
            for (size_t i = 0; i < n_inputs; ++i) U(i,k) = solver.get_solution_value(U_indices[i][k]);
        }
        sigma = solver.get_solution_value(sigma_index);


        // Write solution to files

        {
            ofstream f(file_name_prefix + "X.txt");
            f << X;
        }
        {
            ofstream f(file_name_prefix + "U.txt");
            f << U;
        }

        cout << "sigma   " << sigma << endl;
        cout << "Delta_sigma   " << solver.get_solution_value("Delta_sigma", {}) << endl;
        cout << "==========================================================" << endl;
    }
}