/*
 * 
 * Independent implementation of 
 *     Successive Convexification for 6-DoF Mars Rocket Powered Landing 
 *     with Free-Final-Time (Michael Szmuk, Behcet Acikmese)
 * 
 * https://arxiv.org/abs/1802.03827
 * 
 */

#include "active_model.hpp"
#include "EcosWrapper.hpp"

#include <iostream>
#include <array>
#include <cmath>
#include <ctime>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using std::array;
using std::cout;
using std::endl;


class DiscretizationODE {
private:
    Model::ControlVector u_t, u_t1;
    double sigma, dt;
    Model& model;

public:

    static constexpr size_t n_V_states = 3 + Model::n_states + 2 * Model::n_inputs;
    using state_type = Eigen::Matrix<double, Model::n_states, n_V_states>;

    DiscretizationODE(
        const Model::ControlVector &u_t, 
        const Model::ControlVector &u_t1, 
        const double &sigma, 
        double dt,
        Model& model
    )
    :u_t(u_t)
    ,u_t1(u_t1)
    ,sigma(sigma)
    ,dt(dt)
    ,model(model) {}

    void operator()(const state_type &V, state_type &dVdt, const double t){

        const Model::StateVector &x = V.col(0);
        const Model::ControlVector u = u_t + t / dt * (u_t1 - u_t);

        const double alpha = t / dt;
        const double beta = 1. - alpha;

        const Model::StateMatrix   A_bar  = sigma * model.state_jacobian(x, u);
        const Model::ControlMatrix B_bar  = sigma * model.control_jacobian(x, u);
        const Model::StateVector   f      =         model.ode(x, u);


        Model::StateMatrix Phi_A_xi = V.block<Model::n_states, Model::n_states>(0, 1);
        Model::StateMatrix Phi_A_xi_inverse = Phi_A_xi.inverse();

        size_t cols = 0;

        dVdt.block<Model::n_states,               1>(0, cols) = sigma * f;                                   cols += 1;
        dVdt.block<Model::n_states, Model::n_states>(0, cols) = A_bar * Phi_A_xi;                            cols += Model::n_states;
        dVdt.block<Model::n_states, Model::n_inputs>(0, cols) = Phi_A_xi_inverse * B_bar * alpha;            cols += Model::n_inputs;
        dVdt.block<Model::n_states, Model::n_inputs>(0, cols) = Phi_A_xi_inverse * B_bar * beta;             cols += Model::n_inputs;
        dVdt.block<Model::n_states,               1>(0, cols) = Phi_A_xi_inverse * f;                        cols += 1;
        dVdt.block<Model::n_states,               1>(0, cols) = Phi_A_xi_inverse * (-A_bar * x - B_bar * u);

    }
};

int main() {
    Model model;

    // trajectory points
    constexpr int K = 20;
    const double dt = 1 / double(K-1);

    const double weight_trust_region_sigma = 1.0;
    const double weight_trust_region_xu = 1e-3;
    const double weight_virtual_control = 1e5;

    const size_t n_states = Model::n_states;
    const size_t n_inputs = Model::n_inputs;

    Eigen::Matrix<double, n_states, K> X;
    Eigen::Matrix<double, n_inputs, K> U;


    // START INITIALIZATION
    cout << "Starting initialization." << endl;
    model.initialize<K>(X, U);
    cout << "Initialization finished." << endl;

    // START SUCCESSIVE CONVEXIFICATION
    
    double sigma = model.total_time_guess();


    array<Model::StateMatrix,   (K-1)> A_bar;
    array<Model::ControlMatrix, (K-1)> B_bar;
    array<Model::ControlMatrix, (K-1)> C_bar;
    array<Model::StateVector,   (K-1)> Sigma_bar;
    array<Model::StateVector,   (K-1)> z_bar;


    using namespace boost::numeric::odeint;
    runge_kutta_dopri5<DiscretizationODE::state_type, double, DiscretizationODE::state_type, double, vector_space_algebra> stepper;

    /** Solver setup **/
    EcosWrapper solver;
    {
        solver.create_tensor_variable("X", {n_states, K});
        solver.create_tensor_variable("U", {n_inputs, K});
        solver.create_tensor_variable("nu", {n_states, K-1});
        solver.create_tensor_variable("sigma", {});

        // shortcuts to access solver variables and create parameters
        auto var = [&](const string &name, const vector<size_t> &indices){ return solver.get_variable(name,indices); };
        auto param = [](double &param_value){ return optimization_problem::Parameter(&param_value); };

        // Build linearized model equality constraint
        //   x(k+1) == A x(k) + B u(k) + C u(k+1) + Sigma sigma + z + nu
        // -I x(k+1) + A x(k) + B u(k) + C u(k+1) + Sigma sigma + z + nu == 0

        for (size_t k = 0; k < K-1; k++) {
            for (size_t row_index = 0; row_index < n_states; ++row_index) {

                // -I * x(k+1)
                optimization_problem::AffineExpression eq = (-1.0) * var("X", {row_index, k+1});

                // A * x(k)
                for (size_t col_index = 0; col_index < n_states; ++col_index)
                    eq = eq + param(A_bar.at(k)(row_index, col_index)) * var("X", {col_index, k});

                // B * u(k)
                for (size_t col_index = 0; col_index < n_inputs; ++col_index)
                    eq = eq + param(B_bar.at(k)(row_index, col_index)) * var("U", {col_index, k});

                // C * u(k+1)
                for (size_t col_index = 0; col_index < n_inputs; ++col_index)
                    eq = eq + param(C_bar.at(k)(row_index, col_index)) * var("U", {col_index, k+1});

                // Sigma sigma
                eq = eq + param(Sigma_bar.at(k)(row_index, 0)) * var("sigma", {});

                // z
                eq = eq + param(z_bar.at(k)(row_index, 0));

                // nu
                eq = eq + (1.0) * var("nu", {row_index, 0});

                solver.add_constraint( eq == 0.0 );
            }
        }
    }


    return 0;


    const int iterations = 1;
    for(int it = 1; it < iterations + 1; it++) {
        cout << "Iteration " << it << endl;
        cout << "Calculating new transition matrices." << endl;

        const clock_t begin_time = clock();

        for (int k = 0; k < K-1; k++) {
            DiscretizationODE::state_type V;
            V.setZero();
            V.col(0) = X.col(k);
            V.block<Model::n_states,Model::n_states>(0, 1).setIdentity();

            DiscretizationODE discretizationODE(U.col(k), U.col(k+1), sigma, dt, model);
            integrate_adaptive(make_controlled(1E-12 , 1E-12 , stepper), discretizationODE, V, 0., dt, dt/10.);

            size_t cols = 1;
            A_bar[k]      =            V.block<Model::n_states,Model::n_states>(0, cols);   cols += Model::n_states;
            B_bar[k]      = A_bar[k] * V.block<Model::n_states,Model::n_inputs>(0, cols);   cols += Model::n_inputs;
            C_bar[k]      = A_bar[k] * V.block<Model::n_states,Model::n_inputs>(0, cols);   cols += Model::n_inputs;
            Sigma_bar[k]  = A_bar[k] * V.block<Model::n_states,1>(0, cols);                 cols += 1;
            z_bar[k]      = A_bar[k] * V.block<Model::n_states,1>(0, cols);

            cout << "A_bar " << k << endl;
            cout << A_bar[k] << endl << endl;
            cout << "B_bar " << k << endl;
            cout << B_bar[k] << endl << endl;
            cout << "C_bar " << k << endl;
            cout << C_bar[k] << endl << endl;
            cout << "Sigma_bar " << k << endl;
            cout << Sigma_bar[k] << endl << endl;
            cout << "z_bar " << k << endl;
            cout << z_bar[k] << endl << endl;
        }
        cout << "Transition matrices calculated in " << double( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds." << endl;

        // TODO: Solve problem.

    }
}