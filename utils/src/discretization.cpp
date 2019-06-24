#include "discretization.hpp"

#include <boost/numeric/odeint.hpp>

#include "eigenIntegration.hpp"
#include <eigen3/unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h>

void eulerLinearDiscretization(Model &model,
                               double ts,
                               const Model::state_vector_t &x_eq,
                               const Model::input_vector_t &u_eq,
                               Model::state_matrix_t &A,
                               Model::control_matrix_t &B)
{
    Model::state_matrix_t A_c;
    Model::control_matrix_t B_c;
    model.computeJacobians(x_eq, u_eq, A_c, B_c);

    A = Model::state_matrix_t::Identity() + ts * A_c;
    B = ts * B_c;
}

class ODEExactLinear
{
private:
    Model::state_matrix_t A_c;
    Model::control_matrix_t B_c;
    double ts;

public:
    ODEExactLinear(
        Model::state_matrix_t &A_c,
        Model::control_matrix_t &B_c,
        double ts) : A_c(A_c), B_c(B_c), ts(ts)
    {
    }

    void operator()(const Model::control_matrix_t &B, Model::control_matrix_t &Bdt, const double t) const
    {
        Bdt = (A_c * (ts - t)).exp() * B_c;
    }
};

void exactLinearDiscretization(Model &model,
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
    model.computeJacobians(x_eq, u_eq, A_c, B_c);
    model.computef(x_eq, u_eq, f);

    using namespace boost::numeric::odeint;
    runge_kutta4<Model::control_matrix_t, double, Model::control_matrix_t, double, vector_space_algebra> stepper;
    ODEExactLinear odeExactLinear(A_c, B_c, ts);

    A = (ts * A_c).exp();
    B.setZero();
    integrate_adaptive(stepper, odeExactLinear, B, 0., ts, ts / 10.);
    z = f - A * x_eq - B * u_eq;
}

class ODEMultipleShootingVariableTime
{
private:
    Model::input_vector_t u_t0, u_t1;
    double T, dt;
    Model &model;

public:
    using ode_matrix_t = Eigen::Matrix<double, Model::state_dim, 1 + Model::state_dim + 2 * Model::input_dim + 2>;

    ODEMultipleShootingVariableTime(
        const Model::input_vector_t &u_t0,
        const Model::input_vector_t &u_t1,
        const double &T,
        double dt,
        Model &model)
        : u_t0(u_t0), u_t1(u_t1), T(T), dt(dt), model(model) {}

    void operator()(const ode_matrix_t &V, ode_matrix_t &dVdt, const double t)
    {
        const Model::state_vector_t &x = V.col(0);
        const Model::input_vector_t u = u_t0 + t / dt * (u_t1 - u_t0);

        Model::state_vector_t f;
        Model::state_matrix_t A_bar;
        Model::control_matrix_t B_bar;
        model.computef(x, u, f);
        model.computeJacobians(x, u, A_bar, B_bar);
        A_bar *= T;
        B_bar *= T;

        const Model::state_matrix_t Phi_A_xi = V.block<Model::state_dim, Model::state_dim>(0, 1);
        const Model::state_matrix_t Phi_A_xi_inverse = Phi_A_xi.inverse();

        size_t cols = 0;

        // state
        dVdt.block<Model::state_dim, 1>(0, cols) = T * f;
        cols += 1;

        // A_bar
        dVdt.block<Model::state_dim, Model::state_dim>(0, cols) = A_bar * Phi_A_xi;
        cols += Model::state_dim;

        // B_bar
        const double alpha = (dt - t) / dt;
        dVdt.block<Model::state_dim, Model::input_dim>(0, cols) = Phi_A_xi_inverse * B_bar * alpha;
        cols += Model::input_dim;

        // C_bar
        const double beta = t / dt;
        dVdt.block<Model::state_dim, Model::input_dim>(0, cols) = Phi_A_xi_inverse * B_bar * beta;
        cols += Model::input_dim;

        // S_bar
        dVdt.block<Model::state_dim, 1>(0, cols) = Phi_A_xi_inverse * f;
        cols += 1;

        // z_bar
        dVdt.block<Model::state_dim, 1>(0, cols) = Phi_A_xi_inverse * (-A_bar * x - B_bar * u);

        assert(cols == dVdt.cols() - 1);
    }
};

void multipleShootingVariableTime(
    Model &model,
    double T,
    const Eigen::MatrixXd &X,
    const Eigen::MatrixXd &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &S_bar,
    Model::state_vector_v_t &z_bar)
{
    const size_t K = X.cols();

    const double dt = 1. / double(K - 1);
    using namespace boost::numeric::odeint;
    runge_kutta4<ODEMultipleShootingVariableTime::ode_matrix_t, double, ODEMultipleShootingVariableTime::ode_matrix_t, double, vector_space_algebra> stepper;

    for (size_t k = 0; k < K - 1; k++)
    {
        ODEMultipleShootingVariableTime::ode_matrix_t V;
        V.setZero();
        V.col(0) = X.col(k);
        V.block<Model::state_dim, Model::state_dim>(0, 1).setIdentity();

        ODEMultipleShootingVariableTime odeMultipleShooting(U.col(k), U.col(k + 1), T, dt, model);

        integrate_adaptive(stepper, odeMultipleShooting, V, 0., dt, dt / 4.);

        size_t cols = 1;

        A_bar[k] = V.block<Model::state_dim, Model::state_dim>(0, cols);
        cols += Model::state_dim;

        B_bar[k] = A_bar[k] * V.block<Model::state_dim, Model::input_dim>(0, cols);
        cols += Model::input_dim;

        C_bar[k] = A_bar[k] * V.block<Model::state_dim, Model::input_dim>(0, cols);
        cols += Model::input_dim;

        S_bar[k] = A_bar[k] * V.block<Model::state_dim, 1>(0, cols);
        cols += 1;

        z_bar[k] = A_bar[k] * V.block<Model::state_dim, 1>(0, cols);
    }
}

class ODEMultipleShooting
{
private:
    Model::input_vector_t u_t0, u_t1;
    double dt;
    Model &model;

public:
    using ode_matrix_t = Eigen::Matrix<double, Model::state_dim, 1 + Model::state_dim + 2 * Model::input_dim + 1>;

    ODEMultipleShooting(
        const Model::input_vector_t &u_t0,
        const Model::input_vector_t &u_t1,
        double dt,
        Model &model)
        : u_t0(u_t0), u_t1(u_t1), dt(dt), model(model) {}

    void operator()(const ode_matrix_t &V, ode_matrix_t &dVdt, const double t)
    {
        const Model::state_vector_t &x = V.col(0);
        const Model::input_vector_t u = u_t0 + t / dt * (u_t1 - u_t0);

        Model::state_vector_t f;
        Model::state_matrix_t A_bar;
        Model::control_matrix_t B_bar;
        model.computef(x, u, f);
        model.computeJacobians(x, u, A_bar, B_bar);

        const Model::state_matrix_t Phi_A_xi = V.block<Model::state_dim, Model::state_dim>(0, 1);
        const Model::state_matrix_t Phi_A_xi_inverse = Phi_A_xi.inverse();

        size_t cols = 0;

        // state
        dVdt.block<Model::state_dim, 1>(0, cols) = f;
        cols += 1;

        // A_bar
        dVdt.block<Model::state_dim, Model::state_dim>(0, cols) = A_bar * Phi_A_xi;
        cols += Model::state_dim;

        // B_bar
        const double alpha = (dt - t) / dt;
        dVdt.block<Model::state_dim, Model::input_dim>(0, cols) = Phi_A_xi_inverse * B_bar * alpha;
        cols += Model::input_dim;

        // C_bar
        const double beta = t / dt;
        dVdt.block<Model::state_dim, Model::input_dim>(0, cols) = Phi_A_xi_inverse * B_bar * beta;
        cols += Model::input_dim;

        // z_bar
        dVdt.block<Model::state_dim, 1>(0, cols) = Phi_A_xi_inverse * (f - A_bar * x - B_bar * u);

        assert(cols == dVdt.cols() - 1);
    }
};

void multipleShooting(
    Model &model,
    double T,
    const Eigen::MatrixXd &X,
    const Eigen::MatrixXd &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &z_bar)
{
    const size_t K = X.cols();

    const double dt = T / double(K - 1);
    using namespace boost::numeric::odeint;
    runge_kutta4<ODEMultipleShooting::ode_matrix_t, double, ODEMultipleShooting::ode_matrix_t, double, vector_space_algebra> stepper;

    for (size_t k = 0; k < K - 1; k++)
    {
        ODEMultipleShooting::ode_matrix_t V;
        V.setZero();
        V.col(0) = X.col(k);
        V.block<Model::state_dim, Model::state_dim>(0, 1).setIdentity();

        ODEMultipleShooting odeMultipleShooting(U.col(k), U.col(k + 1), dt, model);

        integrate_adaptive(stepper, odeMultipleShooting, V, 0., dt, dt / 4.);

        size_t cols = 1;

        A_bar[k] = V.block<Model::state_dim, Model::state_dim>(0, cols);
        cols += Model::state_dim;

        B_bar[k] = A_bar[k] * V.block<Model::state_dim, Model::input_dim>(0, cols);
        cols += Model::input_dim;

        C_bar[k] = A_bar[k] * V.block<Model::state_dim, Model::input_dim>(0, cols);
        cols += Model::input_dim;

        z_bar[k] = A_bar[k] * V.block<Model::state_dim, 1>(0, cols);
    }
}