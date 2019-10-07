#include <simulation.hpp>

#include <boost/numeric/odeint.hpp>

#include "eigenIntegration.hpp"

namespace sim
{

ODE::ODE(Model::ptr_t model, double dt, const Model::input_vector_t &u0, const Model::input_vector_t &u1) : model(model), u0(u0), u1(u1), dt(dt) {}

void ODE::operator()(const Model::state_vector_t &f, Model::state_vector_t &dfdt, const double t)
{
    Model::input_vector_t u = u0 + t / dt * (u1 - u0);
    model->computef(f, u, dfdt);
}

void simulate(Model::ptr_t model, double dt,
              const Model::state_vector_t &x0,
              const Model::input_vector_t &u0,
              const Model::input_vector_t &u1,
              Model::state_vector_t &x1)
{
    using namespace boost::numeric::odeint;
    runge_kutta4<Model::state_vector_t, double, Model::state_vector_t, double, vector_space_algebra> stepper;

    ODE ode(model, dt, u0, u1);

    x1 = x0;

    integrate_adaptive(stepper, ode, x1, 0., dt, dt / 5.);
}

} // namespace sim