#include "model_landing_6dof.h"

#include <iostream>
#include <array>
#include <cmath>
#include <ctime>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using std::array;
using std::cout;
using std::endl;

using Model = model_landing_6dof;

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
    constexpr int K = 50;
    const double dt = 1 / double(K-1);

    const size_t n_states = Model::n_states;
    MatrixXd X(n_states, K);

    const size_t n_inputs = Model::n_inputs;
    MatrixXd U(n_inputs, K);


    // START INITIALIZATION
    cout << "Starting initialization." << endl;
    model.initialize(K, X, U);
    cout << "Initialization finished." << endl;

    // START SUCCESSIVE CONVEXIFICATION
    
    double sigma = model.total_time_guess();


    array<Model::StateMatrix,   K> A_bar;
    array<Model::ControlMatrix, K> B_bar;
    array<Model::ControlMatrix, K> C_bar;
    array<Model::StateVector,   K> Sigma_bar;
    array<Model::StateVector,   K> z_bar;


    using namespace boost::numeric::odeint;
    runge_kutta_dopri5<DiscretizationODE::state_type, double, DiscretizationODE::state_type, double, vector_space_algebra> stepper;

    const int iterations = 15;
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

        }
        cout << "Transition matrices calculated in " << double( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds." << endl;

        // TODO: Solve problem.


    }

}