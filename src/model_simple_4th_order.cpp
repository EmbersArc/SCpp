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


void model_simple_4th_order::add_application_constraints(EcosWrapper &solver, int K) {

    auto var = [&](const string &name, const vector<size_t> &indices){ return solver.get_variable(name,indices); };
    auto param = [](double &param_value){ return optimization_problem::Parameter(&param_value); };


    solver.add_constraint( 1.0 * var("X", {0, 0}) + (-1.0) == 0.0 );

    for (size_t i = 1; i < n_states; ++i) {
        solver.add_constraint( 1.0 * var("X", {i, 0}) == 0.0 );
    }

    for (size_t i = 0; i < n_states; ++i) {
        solver.add_constraint( 1.0 * var("X", {i, K-1}) == 0.0 );
    }

    for (size_t k = 0; k < K; ++k) {
        solver.add_constraint( ( 1.0) * var("U", {0, k}) + (1.0) >= (0.0) );
        solver.add_constraint( (-1.0) * var("U", {0, k}) + (1.0) >= (0.0) );
    }

}