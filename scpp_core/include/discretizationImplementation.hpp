#pragma once

#include <boost/numeric/odeint.hpp>

#include "eigenIntegration.hpp"

#include "activeModel.hpp"

namespace scpp::discretization
{

template <bool INTERPOLATE_INPUT, bool VARIABLE_TIME>
class ODE
{
private:
    Model::input_vector_t u_t0, u_t1;
    double time;
    double dt;
    Model::ptr_t model;

public:
    using ode_matrix_t = typename Eigen::Matrix<double, Model::state_dim,
                                                1 + Model::state_dim + Model::input_dim +
                                                    INTERPOLATE_INPUT * Model::input_dim +
                                                    VARIABLE_TIME + 1>;

    ODE(const Model::input_vector_t &u_t0,
        const Model::input_vector_t &u_t1,
        const double &time,
        double dt,
        Model::ptr_t model)
        : u_t0(u_t0), u_t1(u_t1), time(time), dt(dt), model(model) {}

    void operator()(const ode_matrix_t &V, ode_matrix_t &dVdt, const double t);
};

template <bool INTERPOLATE_INPUT, bool VARIABLE_TIME>
void ODE<INTERPOLATE_INPUT, VARIABLE_TIME>::operator()(const ode_matrix_t &V, ode_matrix_t &dVdt, const double t)
{
    const Model::state_vector_t x = V.col(0);

    Model::input_vector_t u;
    if constexpr (INTERPOLATE_INPUT)
    {
        u = u_t0 + t / dt * (u_t1 - u_t0);
    }
    else
    {
        u = u_t0;
    }

    Model::state_vector_t f;
    Model::state_matrix_t A;
    Model::control_matrix_t B;
    model->computef(x, u, f);
    model->computeJacobians(x, u, A, B);

    if constexpr (VARIABLE_TIME)
    {
        A *= time;
        B *= time;
    }

    const Model::state_matrix_t Phi_A_xi = V.template block<Model::state_dim, Model::state_dim>(0, 1);
    const Model::state_matrix_t Phi_A_xi_inverse = Phi_A_xi.inverse();

    size_t cols = 0;

    // state
    if constexpr (VARIABLE_TIME)
    {
        dVdt.template block<Model::state_dim, 1>(0, cols) = time * f;
    }
    else
    {
        dVdt.template block<Model::state_dim, 1>(0, cols) = f;
    }
    cols += 1;

    // A
    dVdt.template block<Model::state_dim, Model::state_dim>(0, cols).noalias() = A * Phi_A_xi;
    cols += Model::state_dim;

    if constexpr (INTERPOLATE_INPUT)
    {
        // B
        const double alpha = (dt - t) / dt;
        dVdt.template block<Model::state_dim, Model::input_dim>(0, cols).noalias() = Phi_A_xi_inverse * B * alpha;
        cols += Model::input_dim;

        // C
        const double beta = t / dt;
        dVdt.template block<Model::state_dim, Model::input_dim>(0, cols).noalias() = Phi_A_xi_inverse * B * beta;
        cols += Model::input_dim;
    }
    else
    {
        // B
        dVdt.template block<Model::state_dim, Model::input_dim>(0, cols).noalias() = Phi_A_xi_inverse * B;
        cols += Model::input_dim;
    }

    if constexpr (VARIABLE_TIME)
    {
        // s
        dVdt.template block<Model::state_dim, 1>(0, cols).noalias() = Phi_A_xi_inverse * f;
        cols += 1;
        // z
        dVdt.template block<Model::state_dim, 1>(0, cols).noalias() = Phi_A_xi_inverse * (-A * x - B * u);
        cols += 1;
    }
    else
    {
        // z
        dVdt.template block<Model::state_dim, 1>(0, cols).noalias() = Phi_A_xi_inverse * (f - A * x - B * u);
        cols += 1;
    }

    assert(cols == ode_matrix_t::ColsAtCompileTime);
}

template <bool INTERPOLATE_INPUT, bool VARIABLE_TIME>
void multipleShootingImplementation(
    Model::ptr_t model,
    trajectory_data_t &td,
    discretization_data_t &dd)
{
    const size_t K = td.n_X();

    using ODEFun = ODE<INTERPOLATE_INPUT, VARIABLE_TIME>;
    using ode_matrix_t = typename ODEFun::ode_matrix_t;

    double dt = 1. / double(K - 1);

    if constexpr (not VARIABLE_TIME)
    {
        dt *= td.t;
    }

    using namespace boost::numeric::odeint;
    runge_kutta_fehlberg78<ode_matrix_t, double, ode_matrix_t, double, vector_space_algebra> stepper;

    for (size_t k = 0; k < K - 1; k++)
    {
        ode_matrix_t V;
        V.col(0) = td.X.at(k);
        V.template block<Model::state_dim, Model::state_dim>(0, 1).setIdentity();
        V.template rightCols<ode_matrix_t::ColsAtCompileTime - 1 - Model::state_dim>().setZero();

        const Model::input_vector_t u0 = td.U[k];
        const Model::input_vector_t u1 = INTERPOLATE_INPUT ? td.U[k + 1] : u0;
        ODEFun odeMultipleShooting(u0, u1, td.t, dt, model);

        integrate_adaptive(stepper, odeMultipleShooting, V, 0., dt, dt / 5.);

        size_t cols = 1;

        dd.A[k] = V.template block<Model::state_dim, Model::state_dim>(0, cols);
        cols += Model::state_dim;

        dd.B[k].noalias() = dd.A[k] * V.template block<Model::state_dim, Model::input_dim>(0, cols);
        cols += Model::input_dim;

        if constexpr (INTERPOLATE_INPUT)
        {
            dd.C[k].noalias() = dd.A[k] * V.template block<Model::state_dim, Model::input_dim>(0, cols);
            cols += Model::input_dim;
        }

        if constexpr (VARIABLE_TIME)
        {
            dd.s[k].noalias() = dd.A[k] * V.template block<Model::state_dim, 1>(0, cols);
            cols += 1;
        }

        dd.z[k].noalias() = dd.A[k] * V.template block<Model::state_dim, 1>(0, cols);
        cols += 1;

        assert(cols == ode_matrix_t::ColsAtCompileTime);
    }
}

} // namespace scpp::discretization