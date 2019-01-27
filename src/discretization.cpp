
#include "discretization.hpp"

class DiscretizationODE
{
  private:
    Model::input_vector_t u_t0, u_t1;
    double sigma, dt;
    Model &model;

  public:
    using state_type = Eigen::Matrix<double, Model::state_dim_, 1 + Model::state_dim_ + 2 * Model::input_dim_ + 2>;

    DiscretizationODE(
        const Model::input_vector_t &u_t0,
        const Model::input_vector_t &u_t1,
        const double &sigma,
        double dt,
        Model &model)
        : u_t0(u_t0), u_t1(u_t1), sigma(sigma), dt(dt), model(model) {}

    void operator()(const state_type &V, state_type &dVdt, const double t)
    {
        const Model::state_vector_t &x = V.col(0);
        const Model::input_vector_t u = u_t0 + t / dt * (u_t1 - u_t0);

        Model::state_vector_t f;
        Model::state_matrix_t A_bar;
        Model::control_matrix_t B_bar;
        model.computef(x, u, f);
        model.computeJacobians(x, u, A_bar, B_bar);
        A_bar *= sigma;
        B_bar *= sigma;

        const Model::state_matrix_t Phi_A_xi = V.block<Model::state_dim_, Model::state_dim_>(0, 1);
        const Model::state_matrix_t Phi_A_xi_inverse = Phi_A_xi.inverse();

        size_t cols = 0;

        dVdt.block<Model::state_dim_, 1>(0, cols) = sigma * f;
        cols += 1;

        dVdt.block<Model::state_dim_, Model::state_dim_>(0, cols) = A_bar * Phi_A_xi;
        cols += Model::state_dim_;

        const double alpha = (dt - t) / dt;
        dVdt.block<Model::state_dim_, Model::input_dim_>(0, cols) = Phi_A_xi_inverse * B_bar * alpha;
        cols += Model::input_dim_;

        const double beta = t / dt;
        dVdt.block<Model::state_dim_, Model::input_dim_>(0, cols) = Phi_A_xi_inverse * B_bar * beta;
        cols += Model::input_dim_;

        dVdt.block<Model::state_dim_, 1>(0, cols) = Phi_A_xi_inverse * f;
        cols += 1;

        dVdt.block<Model::state_dim_, 1>(0, cols) = Phi_A_xi_inverse * (-A_bar * x - B_bar * u);
    }
};

void calculateDiscretization(
    Model &model,
    double &sigma,
    const Eigen::MatrixXd &X,
    const Eigen::MatrixXd &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &Sigma_bar,
    Model::state_vector_v_t &z_bar)
{
    const size_t K = X.cols();

    const double dt = 1. / double(K - 1);
    using namespace boost::numeric::odeint;
    runge_kutta4<DiscretizationODE::state_type, double, DiscretizationODE::state_type, double, vector_space_algebra> stepper;

    for (size_t k = 0; k < K - 1; k++)
    {
        DiscretizationODE::state_type V;
        V.setZero();
        V.col(0) = X.col(k);
        V.block<Model::state_dim_, Model::state_dim_>(0, 1).setIdentity();

        DiscretizationODE discretizationODE(U.col(k), U.col(k + 1), sigma, dt, model);

        integrate_adaptive(stepper, discretizationODE, V, 0., dt, dt / 5.);

        size_t cols = 1;

        A_bar[k] = V.block<Model::state_dim_, Model::state_dim_>(0, cols);
        cols += Model::state_dim_;

        B_bar[k] = A_bar[k] * V.block<Model::state_dim_, Model::input_dim_>(0, cols);
        cols += Model::input_dim_;

        C_bar[k] = A_bar[k] * V.block<Model::state_dim_, Model::input_dim_>(0, cols);
        cols += Model::input_dim_;

        Sigma_bar[k] = A_bar[k] * V.block<Model::state_dim_, 1>(0, cols);
        cols += 1;

        z_bar[k] = A_bar[k] * V.block<Model::state_dim_, 1>(0, cols);
    }
}