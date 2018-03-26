#include "check_jacobian.h"

TestODE::TestODE(){
    // linear interpolation between random input for t = 0 and t = 1
    U0.setRandom().normalize();
    U0 = U0.cwiseAbs() * 0.1;
    U1.setRandom().normalize();
    U1 = U1.cwiseAbs();
}


void TestODE::operator()(const TestStateVector &X, TestStateVector &dXdt, const double t){
    Vector3d U = U0 + t * (U1 - U0);

    // x' = f(x,u)
    dXdt.head(Model::n_states) = model.ode(X.head(Model::n_states), U);
    // x' = A(x,u) * x + B(x,u) * u
    dXdt.tail(Model::n_states) = model.state_jacobian(X.head(Model::n_states), U) * X.tail(Model::n_states)
                                 + model.control_jacobian(X.head(Model::n_states), U) * U ;
}

bool test_model(Model::StateVector V_init) {
    using namespace boost::numeric::odeint;

    TestStateVector V;
    V << V_init, V_init;

    TestODE testODE;

    runge_kutta_dopri5<TestStateVector, double, TestStateVector, double, vector_space_algebra> stepper;
    integrate_adaptive(make_controlled(1E-12, 1E-12, stepper), testODE, V, 0., 1., 1. / 10000.);

    double square_difference = (V.head(Model::n_states) - V.tail(Model::n_states)).squaredNorm();
    std::cout << square_difference << std::endl;

    // true if Jacobian is consistent
    return square_difference < 1e-20;

}

