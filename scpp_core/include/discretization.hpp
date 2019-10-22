#pragma once

#include "discretitation_ode.hpp"

#include <boost/numeric/odeint.hpp>

#include "eigenIntegration.hpp"
#include <eigen3/unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h>

namespace scpp::discretization
{

void exactLinearDiscretization(Model::ptr_t model,
                               double ts,
                               const Model::state_vector_t &x_eq,
                               const Model::input_vector_t &u_eq,
                               Model::state_matrix_t &A,
                               Model::control_matrix_t &B,
                               Model::state_vector_t &z);

template <InputType INPUT_TYPE, TimeType TIME_TYPE>
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

    typedef ODE<INPUT_TYPE, TIME_TYPE> ODEFun;

    double dt = 1. / double(K - 1);

    if constexpr (TIME_TYPE == TimeType::fixed)
    {
        dt *= T;
    }

    using namespace boost::numeric::odeint;
    runge_kutta4<typename ODEFun::ode_matrix_t, double, typename ODEFun::ode_matrix_t, double, vector_space_algebra> stepper;

    for (size_t k = 0; k < K - 1; k++)
    {
        typename ODEFun::ode_matrix_t V;
        V.template setZero();
        V.col(0) = X.at(k);
        V.template block<Model::state_dim, Model::state_dim>(0, 1).setIdentity();

        ODEFun odeMultipleShooting(U[k], U[k + 1], T, dt, model);

        integrate_adaptive(stepper, odeMultipleShooting, V, 0., dt, dt / 3.);

        size_t cols = 1;

        A_bar[k] = V.template block<Model::state_dim, Model::state_dim>(0, cols);
        cols += Model::state_dim;

        B_bar[k].noalias() = A_bar[k] * V.template block<Model::state_dim, Model::input_dim>(0, cols);
        cols += Model::input_dim;

        if constexpr (INPUT_TYPE == InputType::interpolated)
        {
            C_bar[k].noalias() = A_bar[k] * V.template block<Model::state_dim, Model::input_dim>(0, cols);
            cols += Model::input_dim;
        }

        if constexpr (TIME_TYPE == TimeType::variable)
        {
            S_bar[k].noalias() = A_bar[k] * V.template block<Model::state_dim, 1>(0, cols);
            cols += 1;
        }

        z_bar[k].noalias() = A_bar[k] * V.template block<Model::state_dim, 1>(0, cols);
        cols += 1;

        assert(cols == ODEFun::ode_matrix_t::ColsAtCompileTime);
    }
}

} // namespace scpp::discretization