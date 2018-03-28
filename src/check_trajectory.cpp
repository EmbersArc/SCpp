#include "check_trajectory.h"


void ODE::operator()(const Model::StateVector &X, Model::StateVector &dXdt, const double t){
    Model::ControlVector U = U0 + t / dt * (U1 - U0);
    dXdt = model.ode(X,U);
}

void ODE::update_input(const Model::ControlVector &U0, const Model::ControlVector &U1) {
    this->U0 = U0;
    this->U1 = U1;
}

