
#include "EcosWrapper.hpp"
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace optimization_problem;

#define K (50)
#define n_states (14)
#define n_inputs (3)

int main() {

    double weight_trust_region_xu = 1.0;


    Eigen::Matrix<double, n_states, K> X;
    Eigen::Matrix<double, n_inputs, K> U;
    X.setRandom();
    U.setRandom();
    X*=100;


    EcosWrapper solver;


    solver.create_tensor_variable("X", {n_states, K}); // states
    solver.create_tensor_variable("U", {n_inputs, K}); // inputs
    solver.create_tensor_variable("Delta", {K}); // squared change of the stacked [ x(k), u(k) ] vector
    solver.create_tensor_variable("norm2_Delta", {}); // 2-norm of the Delta(k) variables

    // shortcuts to access solver variables and create parameters
    auto var = [&](const string &name, const vector<size_t> &indices){ return solver.get_variable(name,indices); };
    auto param = [](double &param_value){ return optimization_problem::Parameter(&param_value); };
    auto param_fn = [](std::function<double()> callback){ return optimization_problem::Parameter(callback); };


    
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


    solver.compile_problem_structure();
    solver.solve_problem();

    // Read solution
    Eigen::Matrix<double, n_states, K> X2;
    Eigen::Matrix<double, n_inputs, K> U2;
    for (size_t k = 0; k < K; k++) {
        for (size_t i = 0; i < n_states; ++i) X2(i,k) = solver.get_solution_value("X", {i,k});
        for (size_t i = 0; i < n_inputs; ++i) U2(i,k) = solver.get_solution_value("U", {i,k});
    }


    // Should be all zeros
    /*cout << endl << "dX" << endl;
    cout << X-X2 << endl;

    cout << endl << "dU" << endl;
    cout << U-U2 << endl;*/

    cout << endl << "max(abs(dX))" << endl;
    cout << (X-X2).array().abs().maxCoeff() << endl;

    cout << endl << "max(abs(dU))" << endl;
    cout << (U-U2).array().abs().maxCoeff() << endl;
}