#include "check_trajectory.h"


void ODE::operator()(const Model::StateVector &X, Model::StateVector &dXdt, const double t){
    Model::ControlVector U = U0 + t / dt * (U1 - U0);
    dXdt = model.ode(X,U);
}

void ODE::update_input(const Model::ControlVector &U0, const Model::ControlVector &U1) {
    this->U0 = U0;
    this->U1 = U1;
}

void check_trajectory(
        Eigen::Matrix<double, Model::n_states, K> &X_in,
        Eigen::Matrix<double, Model::n_inputs, K> &U,
        double sigma
) {
    using namespace boost::numeric::odeint;

    MatrixXd X_exact(Model::n_states, K);
    X_exact.col(0) = X_in.col(0);

    Model::StateVector x = X_in.col(0);

    runge_kutta_dopri5<Model::StateVector , double, Model::StateVector , double, vector_space_algebra> stepper;


    double dt = sigma / (K - 1);
    ODE ode(dt);

    for(int k = 1; k < K; k++){
        ode.update_input(U.col(k-1), U.col(k));
        integrate_adaptive(make_controlled(1e-12, 1e-12, stepper), ode, x, 0., dt, dt/15);
        X_exact.col(k) = x;
    }

    double total_absolute_error = (X_exact - X_in).cwiseAbs().sum();
    std::cout << std::endl;
    std::cout << "Absolute errors: " << std::endl;
    std::cout << X_exact - X_in << std::endl;
    std::cout << std::endl << "Total absolute error: " << total_absolute_error << std::endl << std::endl;

}