#include "model_simple_4th_order.hpp"


model_simple_4th_order::StateMatrix model_simple_4th_order::state_jacobian(const StateVector &x, const ControlVector &u) {
    StateMatrix A;
    A.setZero();
    A(0, 1) = 1.0;
    A(1, 2) = 1.0;
    A(2, 3) = 1.0;
    return A;
}



model_simple_4th_order::ControlMatrix model_simple_4th_order::control_jacobian(const StateVector &x, const ControlVector &u) {
    ControlMatrix B;
    B.setZero();
    B(3,0) = 1.0;
    return B; 
}



model_simple_4th_order::StateVector model_simple_4th_order::ode(const StateVector &x, const ControlVector &u) {
    StateVector f;
    f << x[1], x[2], x[3], u[0];
    return f;
}
