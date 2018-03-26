#include "check_jacobian.h"

bool test_model() {
    Model model;

    Model::StateVector X, dX;
    Model::ControlVector U, dU;


    double lin, ode;

    double epsilon = 1e-20;

    for(int i = 0; i < Model::n_states; i++) {
        for(int j = 0; j < Model::n_states; j++) {
            X.setOnes();
            dX.setZero();
            dX(j) = epsilon;
            lin = epsilon * model.state_jacobian(X, U)(i, j);
            ode = (model.ode(X + dX, U) - model.ode(X, U))(i);
            std::cout << abs(lin - ode) << std::endl;

            if(abs(lin - ode) > 1e-15){
                return false;
            };
        }
    }

    for(int i = 0; i < Model::n_states; i++) {
        for(int j = 0; j < Model::n_inputs; j++) {
            U.setConstant(1. / sqrt(3));
            dU.setZero();
            dU(j) = epsilon;

            lin = epsilon * model.control_jacobian(X, U)(i, j);
            ode = (model.ode(X, U + dU) - model.ode(X, U))(i);
            std::cout << abs(lin - ode) << std::endl;

            if(abs(lin - ode) > 1e-15){
                return false;
            };
        }
    }

    return true;

}

