#include "check_jacobian.h"
#include "active_model.hpp"

bool check_jacobian() {
    Model model;

    Model::StateVector X, dX;
    Model::ControlVector U, dU;
    X.setRandom();
    U.setRandom().normalize();

    double epsilon = 1e-6;

    auto state_jacobian = model.state_jacobian(X, U);

    for(int i = 0; i < Model::n_states; i++) {
        for(int j = 0; j < Model::n_states; j++) {
            dX.setZero();
            dX(j) = epsilon;
            double analytic_derivative = state_jacobian(i, j);
            double numeric_derivative = (model.ode(X + dX, U) - model.ode(X, U))(i) / epsilon;

            if(abs(analytic_derivative - numeric_derivative) > 1e-5){
                return false;
            };
        }
    }

    auto control_jacobian = model.control_jacobian(X, U);

    for(int i = 0; i < Model::n_states; i++) {
        for(int j = 0; j < Model::n_inputs; j++) {
            dU.setZero();
            dU(j) = epsilon;

            double analytic_derivative = control_jacobian(i, j);
            double numeric_derivative = (model.ode(X, U + dU) - model.ode(X, U))(i) / epsilon;

            if(abs(analytic_derivative - numeric_derivative) > 1e-5){
                return false;
            };
        }
    }

    return true;

}

