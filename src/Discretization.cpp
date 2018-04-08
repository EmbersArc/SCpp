#include "Discretization.hpp"
#include "constants.hpp"

#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>
#include <boost/numeric/odeint.hpp>


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



void calculate_discretization (
    Model &model,
    double &sigma,
    Eigen::Matrix<double, Model::n_states, K> &X,
    Eigen::Matrix<double, Model::n_inputs, K> &U,
    array<Model::StateMatrix,   (K-1)> &A_bar,
    array<Model::ControlMatrix, (K-1)> &B_bar,
    array<Model::ControlMatrix, (K-1)> &C_bar,
    array<Model::StateVector,   (K-1)> &Sigma_bar,
    array<Model::StateVector,   (K-1)> &z_bar
) {

    const double dt = 1 / double(K-1);
    using namespace boost::numeric::odeint;
    runge_kutta4<DiscretizationODE::state_type, double, DiscretizationODE::state_type, double, vector_space_algebra> stepper;


    for (size_t k = 0; k < K-1; k++) {
        DiscretizationODE::state_type V;
        V.setZero();
        V.col(0) = X.col(k);
        V.block<Model::n_states,Model::n_states>(0, 1).setIdentity();

        DiscretizationODE discretizationODE(U.col(k), U.col(k+1), sigma, dt, model);
        integrate_n_steps( stepper , discretizationODE , V , 0. , dt/10.0 , 10 );

        size_t cols = 1;
        A_bar[k]      =            V.block<Model::n_states,Model::n_states>(0, cols);   cols += Model::n_states;
        B_bar[k]      = A_bar[k] * V.block<Model::n_states,Model::n_inputs>(0, cols);   cols += Model::n_inputs;
        C_bar[k]      = A_bar[k] * V.block<Model::n_states,Model::n_inputs>(0, cols);   cols += Model::n_inputs;
        Sigma_bar[k]  = A_bar[k] * V.block<Model::n_states,1>(0, cols);                 cols += 1;
        z_bar[k]      = A_bar[k] * V.block<Model::n_states,1>(0, cols);

        /*cout << "A_bar " << k << endl;
        cout << A_bar[k] << endl << endl;
        cout << "B_bar " << k << endl;
        cout << B_bar[k] << endl << endl;
        cout << "C_bar " << k << endl;
        cout << C_bar[k] << endl << endl;
        cout << "Sigma_bar " << k << endl;
        cout << Sigma_bar[k] << endl << endl;
        cout << "z_bar " << k << endl;
        cout << z_bar[k] << endl << endl;*/
    }


}



