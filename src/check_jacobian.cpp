#include "check_jacobian.h"
#include "active_model.hpp"
#include <cmath>

void check_jacobian(
    // inputs
    double epsilon, double random_radius, 
    // outputs
    double &max_absolute_error, double &max_relative_error
){
    Model model;

    Model::StateVector X, dX;
    Model::ControlVector U, dU;
    X.setRandom();
    U.setRandom();

    X *= random_radius;
    U *= random_radius;

    auto state_jacobian = model.state_jacobian(X, U);
    auto control_jacobian = model.control_jacobian(X, U);

    max_absolute_error = 0;
    max_relative_error = 0;

    for(size_t i = 0; i < Model::n_states; i++) {
        for(size_t j = 0; j < Model::n_states; j++) {
            dX.setZero();
            dX(j) = epsilon;
            double analytic_derivative = state_jacobian(i, j);
            double numeric_derivative = (model.ode(X + dX, U) - model.ode(X, U))(i) / epsilon;

            double absolute_error = fabs(analytic_derivative - numeric_derivative);
            double relative_error = absolute_error / (1e-30 + fmax(fabs(analytic_derivative), fabs(numeric_derivative)));

            max_absolute_error = fmax(max_absolute_error, absolute_error);
            max_relative_error = fmax(max_relative_error, relative_error);
        }
    }

    for(size_t i = 0; i < Model::n_states; i++) {
        for(size_t j = 0; j < Model::n_inputs; j++) {
            dU.setZero();
            dU(j) = epsilon;

            double analytic_derivative = control_jacobian(i, j);
            double numeric_derivative = (model.ode(X, U + dU) - model.ode(X, U))(i) / epsilon;

            double absolute_error = fabs(analytic_derivative - numeric_derivative);
            double relative_error = absolute_error / (1e-30 + fmax(fabs(analytic_derivative), fabs(numeric_derivative)));

            max_absolute_error = fmax(max_absolute_error, absolute_error);
            max_relative_error = fmax(max_relative_error, relative_error);
        }
    }

}

