#include "model_landing_3dof.h"



model_landing_3dof::StateVector model_landing_3dof::ode(const StateVector &x, const ControlVector &u) {

    const double throttle = u(0, 0);
    const double gimbalAngle = u(1, 0);
    const double vx = x(2, 0);
    const double vy = x(3, 0);
    const double theta = x(4, 0);
    const double dtheta = x(5, 0);
    const double x0 = TWR*g*throttle;
    const double x1 = gimbalAngle + theta;

    StateVector f;
    f(0, 0) = vx;
    f(1, 0) = vy;
    f(2, 0) = x0*sin(x1);
    f(3, 0) = g*(TWR*throttle*cos(x1) - 1);
    f(4, 0) = dtheta;
    f(5, 0) = -rTB*x0*sin(gimbalAngle)/(rG*rG);


    return f;
}

model_landing_3dof::StateMatrix model_landing_3dof::state_jacobian(const StateVector &x, const ControlVector &u) {

    const double throttle = u(0, 0);
    const double gimbalAngle = u(1, 0);
    const double theta = x(4, 0);
    const double x0 = TWR*g*throttle;
    const double x1 = gimbalAngle + theta;

    StateMatrix A;
    A.setZero();
    A(0, 2) = 1;
    A(1, 3) = 1;
    A(2, 4) = x0*cos(x1);
    A(3, 4) = -x0*sin(x1);
    A(4, 5) = 1;


    return A;
}



model_landing_3dof::ControlMatrix model_landing_3dof::control_jacobian(const StateVector &x, const ControlVector &u) {

    const double throttle = u(0, 0);
    const double gimbalAngle = u(1, 0);
    const double theta = x(4, 0);
    const double x0 = TWR*g;
    const double x1 = gimbalAngle + theta;
    const double x2 = x0*sin(x1);
    const double x3 = x0*cos(x1);
    const double x4 = TWR*g*rTB/(rG*rG);

    ControlMatrix B;
    B.setZero();
    B(2, 0) = x2;
    B(2, 1) = throttle*x3;
    B(3, 0) = x3;
    B(3, 1) = -throttle*x2;
    B(5, 0) = -x4*sin(gimbalAngle);
    B(5, 1) = -throttle*x4*cos(gimbalAngle);


    return B; 
}



void model_landing_3dof::add_application_constraints(EcosWrapper &solver, size_t K) {
    // TODO
}

model_landing_3dof::StateVector model_landing_3dof::get_random_state() {
    StateVector X;
    X.setRandom();
    X *= 10;
    return X;
}

model_landing_3dof::ControlVector model_landing_3dof::get_random_input() {
    ControlVector U;
    U.setRandom();
    return U;
}