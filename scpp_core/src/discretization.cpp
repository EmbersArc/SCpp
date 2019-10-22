#include "discretization.hpp"
#include "discretitation_ode.hpp"

#include <boost/numeric/odeint.hpp>

#include "eigenIntegration.hpp"
#include <eigen3/unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h>

namespace scpp
{
namespace discretization
{

void multipleShooting(
    Model::ptr_t model,
    double T,
    const Model::state_vector_v_t &X,
    const Model::input_vector_v_t &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &S_bar,
    Model::state_vector_v_t &z_bar)
{
    const size_t K = X.size();

    using ODEMultipleShootingVariableTime = ODE<InputType::interpolated, TimeType::variable>;

    const double dt = 1. / double(K - 1);
    using namespace boost::numeric::odeint;
    runge_kutta4<ODEMultipleShootingVariableTime::ode_matrix_t,
                 double,
                 ODEMultipleShootingVariableTime::ode_matrix_t,
                 double,
                 vector_space_algebra>
        stepper;

    for (size_t k = 0; k < K - 1; k++)
    {
        ODEMultipleShootingVariableTime::ode_matrix_t V;
        V.setZero();
        V.col(0) = X.at(k);
        V.block<Model::state_dim, Model::state_dim>(0, 1).setIdentity();

        ODEMultipleShootingVariableTime odeMultipleShooting(U[k], U[k + 1], T, dt, model);

        integrate_adaptive(stepper, odeMultipleShooting, V, 0., dt, dt / 3.);

        size_t cols = 1;

        A_bar[k] = V.block<Model::state_dim, Model::state_dim>(0, cols);
        cols += Model::state_dim;

        B_bar[k].noalias() = A_bar[k] * V.block<Model::state_dim, Model::input_dim>(0, cols);
        cols += Model::input_dim;

        C_bar[k].noalias() = A_bar[k] * V.block<Model::state_dim, Model::input_dim>(0, cols);
        cols += Model::input_dim;

        S_bar[k].noalias() = A_bar[k] * V.block<Model::state_dim, 1>(0, cols);
        cols += 1;

        z_bar[k].noalias() = A_bar[k] * V.block<Model::state_dim, 1>(0, cols);
    }
}

void multipleShooting(
    Model::ptr_t model,
    double T,
    const Model::state_vector_v_t &X,
    const Model::input_vector_v_t &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &z_bar)
{
    const size_t K = X.size();

    using ODEMultipleShooting = ODE<InputType::interpolated, TimeType::fixed>;

    const double dt = T / double(K - 1);
    using namespace boost::numeric::odeint;
    runge_kutta4<ODEMultipleShooting::ode_matrix_t,
                 double,
                 ODEMultipleShooting::ode_matrix_t,
                 double,
                 vector_space_algebra>
        stepper;

    for (size_t k = 0; k < K - 1; k++)
    {
        ODEMultipleShooting::ode_matrix_t V;
        V.setZero();
        V.col(0) = X.at(k);
        V.block<Model::state_dim, Model::state_dim>(0, 1).setIdentity();

        ODEMultipleShooting odeMultipleShooting(U[k], U[k + 1], dt, model);

        integrate_adaptive(stepper, odeMultipleShooting, V, 0., dt, dt / 3.);

        size_t cols = 1;

        A_bar[k] = V.block<Model::state_dim, Model::state_dim>(0, cols);
        cols += Model::state_dim;

        B_bar[k].noalias() = A_bar[k] * V.block<Model::state_dim, Model::input_dim>(0, cols);
        cols += Model::input_dim;

        C_bar[k].noalias() = A_bar[k] * V.block<Model::state_dim, Model::input_dim>(0, cols);
        cols += Model::input_dim;

        z_bar[k].noalias() = A_bar[k] * V.block<Model::state_dim, 1>(0, cols);
    }
}

void exactLinearDiscretization(Model::ptr_t model,
                               double ts,
                               const Model::state_vector_t &x_eq,
                               const Model::input_vector_t &u_eq,
                               Model::state_matrix_t &A,
                               Model::control_matrix_t &B,
                               Model::state_vector_t &z)
{
    Model::state_matrix_t A_c;
    Model::control_matrix_t B_c;
    Model::state_vector_t f;
    model->computeJacobians(x_eq, u_eq, A_c, B_c);
    model->computef(x_eq, u_eq, f);

    Eigen::MatrixXd E;
    E.resize(Model::state_dim + Model::input_dim, Model::state_dim + Model::input_dim);
    E.setZero();
    E.topLeftCorner<Model::state_dim, Model::state_dim>() = A_c;
    E.topRightCorner<Model::state_dim, Model::input_dim>() = B_c;
    Eigen::MatrixXd expE = (E * ts).exp();

    A = expE.topLeftCorner<Model::state_dim, Model::state_dim>();
    B = expE.topRightCorner<Model::state_dim, Model::input_dim>();

    E.resize(Model::state_dim + 1, Model::state_dim + 1);
    E.setZero();
    E.topLeftCorner<Model::state_dim, Model::state_dim>() = A_c;
    E.topRightCorner<Model::state_dim, 1>() = f - A_c * x_eq - B_c * u_eq;
    expE = (E * ts).exp();

    z = expE.topRightCorner<Model::state_dim, 1>();
}

} // namespace discretization
} // namespace scpp