#pragma once

#include "activeModel.hpp"

namespace scpp
{
namespace discretization
{

enum InputType : size_t
{
    constant = 0,
    interpolated = 1
};

enum TimeType : size_t
{
    fixed = 0,
    variable = 1
};

template <InputType INPUT_TYPE, TimeType TIME_TYPE>
class ODE
{
private:
    Model::input_vector_t u_t0, u_t1;
    double T;
    double dt;
    Model::ptr_t model;

public:
    using ode_matrix_t = Eigen::Matrix<double, Model::state_dim,
                                       1 + Model::state_dim + Model::input_dim +
                                           INPUT_TYPE * Model::input_dim +
                                           TIME_TYPE + 1>;

    ODE(const Model::input_vector_t &u_t0, double dt, Model::ptr_t model)
        : u_t0(u_t0), u_t1(0), T(0), dt(dt), model(model)
    {
        static_assert(INPUT_TYPE == InputType::constant and TIME_TYPE == TimeType::fixed);
    }

    ODE(const Model::input_vector_t &u_t0,
        const double &T,
        double dt,
        Model::ptr_t model)
        : u_t0(u_t0), u_t1(0), T(T), dt(dt), model(model)
    {
        static_assert(INPUT_TYPE == InputType::constant and TIME_TYPE == TimeType::variable);
    }

    ODE(const Model::input_vector_t &u_t0,
        const Model::input_vector_t &u_t1,
        double dt,
        Model::ptr_t model)
        : u_t0(u_t0), u_t1(u_t1), T(0), dt(dt), model(model)
    {
        static_assert(INPUT_TYPE == InputType::fixed and TIME_TYPE == TimeType::fixed);
    }

    ODE(const Model::input_vector_t &u_t0,
        const Model::input_vector_t &u_t1,
        const double &T,
        double dt,
        Model::ptr_t model)
        : u_t0(u_t0), u_t1(u_t1), T(T), dt(dt), model(model)
    {
        static_assert(INPUT_TYPE == InputType::interpolated and TIME_TYPE == TimeType::variable);
    }

    void operator()(const ode_matrix_t &V, ode_matrix_t &dVdt, const double t);
};

template <InputType INPUT_TYPE, TimeType TIME_TYPE>
void ODE<INPUT_TYPE, TIME_TYPE>::operator()(const ode_matrix_t &V, ode_matrix_t &dVdt, const double t)
{
    const Model::state_vector_t x = V.col(0);

    Model::input_vector_t u;
    if constexpr (INPUT_TYPE == InputType::interpolated)
    {
        u = u_t0 + t / dt * (u_t1 - u_t0);
    }
    else
    {
        u = u_t0;
    }

    Model::state_vector_t f;
    Model::state_matrix_t A_bar;
    Model::control_matrix_t B_bar;
    model->computef(x, u, f);
    model->computeJacobians(x, u, A_bar, B_bar);

    if constexpr (TIME_TYPE == TimeType::variable)
    {
        A_bar *= T;
        B_bar *= T;
    }

    const Model::state_matrix_t Phi_A_xi = V.template block<Model::state_dim, Model::state_dim>(0, 1);
    const Model::state_matrix_t Phi_A_xi_inverse = Phi_A_xi.inverse();

    size_t cols = 0;

    // state
    if constexpr (TIME_TYPE == TimeType::fixed)
    {
        dVdt.template block<Model::state_dim, 1>(0, cols) = f;
    }
    else
    {
        dVdt.template block<Model::state_dim, 1>(0, cols) = T * f;
    }
    cols += 1;

    // A_bar
    dVdt.template block<Model::state_dim, Model::state_dim>(0, cols).noalias() = A_bar * Phi_A_xi;
    cols += Model::state_dim;

    if constexpr (INPUT_TYPE == InputType::constant)
    {
        // B_bar
        dVdt.template block<Model::state_dim, Model::input_dim>(0, cols).noalias() = Phi_A_xi_inverse * B_bar;
        cols += Model::input_dim;
    }
    else
    {
        // B_bar
        const double alpha = (dt - t) / dt;
        dVdt.template block<Model::state_dim, Model::input_dim>(0, cols).noalias() = Phi_A_xi_inverse * B_bar * alpha;
        cols += Model::input_dim;

        // C_bar
        const double beta = t / dt;
        dVdt.template block<Model::state_dim, Model::input_dim>(0, cols).noalias() = Phi_A_xi_inverse * B_bar * beta;
        cols += Model::input_dim;
    }

    if constexpr (TIME_TYPE == TimeType::fixed)
    {
        // z_bar
        dVdt.template block<Model::state_dim, 1>(0, cols).noalias() = Phi_A_xi_inverse * (f - A_bar * x - B_bar * u);
        cols += 1;
    }
    else
    {
        // S_bar
        dVdt.template block<Model::state_dim, 1>(0, cols).noalias() = Phi_A_xi_inverse * f;
        cols += 1;
        // z_bar
        dVdt.template block<Model::state_dim, 1>(0, cols).noalias() = Phi_A_xi_inverse * (-A_bar * x - B_bar * u);
        cols += 1;
    }

    assert(cols == size_t(dVdt.cols()));
}

} // namespace discretization
} // namespace scpp