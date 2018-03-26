#include "check_jacobian.h"

bool test_model() {
    Model model;

    Model::StateVector X, dX;
    Model::ControlVector U, dU;
    X.setRandom();
    U.setRandom().normalize();


    double lin, ode;

    double epsilon = 1e-10;

    for(int i = 0; i < Model::n_states; i++) {
        for(int j = 0; j < Model::n_states; j++) {
            dX.setZero();
            dX(j) = epsilon;
            lin = model.state_jacobian(X, U)(i, j);
            ode = (model.ode(X + dX, U) - model.ode(X, U))(i) / epsilon;
            std::cout << abs(lin - ode) << std::endl;

            if(abs(lin - ode) > 1e-5){
                return false;
            };
        }
    }

    for(int i = 0; i < Model::n_states; i++) {
        for(int j = 0; j < Model::n_inputs; j++) {
            dU.setZero();
            dU(j) = epsilon;

            lin = model.control_jacobian(X, U)(i, j);
            ode = (model.ode(X, U + dU) - model.ode(X, U))(i) / epsilon;
            std::cout << abs(lin - ode) << std::endl;

            if(abs(lin - ode) > 1e-5){
                return false;
            };
        }
    }

    return true;

}

