#include "discretization.hpp"
#include "constants.hpp"

class DiscretizationODE
{
  private:
    Model::input_vector_t u_t, u_t1;
    double sigma, dt;
    Model &model;

  public:
    static constexpr size_t n_V_states = 3 + Model::state_dim_ + 2 * Model::input_dim_;
    using state_type = Eigen::Matrix<double, Model::state_dim_, n_V_states>;

    DiscretizationODE(
        const Model::input_vector_t &u_t,
        const Model::input_vector_t &u_t1,
        const double &sigma,
        double dt,
        Model &model)
        : u_t(u_t), u_t1(u_t1), sigma(sigma), dt(dt), model(model) {}

    void operator()(const state_type &V, state_type &dVdt, const double t)
    {
        const Model::state_vector_t x = V.col(0);
        const Model::input_vector_t u = u_t + t / dt * (u_t1 - u_t);

        Model::state_vector_t f;
        model.computef(x, u, f);
        Model::state_matrix_t A_bar;
        model.computeA(x, u, A_bar);
        A_bar *= sigma;
        Model::control_matrix_t B_bar;
        model.computeB(x, u, B_bar);
        B_bar *= sigma;

        Model::state_matrix_t Phi_A_xi = V.block<Model::state_dim_, Model::state_dim_>(0, 1);
        Model::state_matrix_t Phi_A_xi_inverse = Phi_A_xi.inverse();

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

void calculate_discretization(
    Model &model,
    double &sigma,
    Eigen::Matrix<double, Model::state_dim_, K> &X,
    Eigen::Matrix<double, Model::input_dim_, K> &U,
    array<Model::state_matrix_t, (K - 1)> &A_bar,
    array<Model::control_matrix_t, (K - 1)> &B_bar,
    array<Model::control_matrix_t, (K - 1)> &C_bar,
    array<Model::state_vector_t, (K - 1)> &Sigma_bar,
    array<Model::state_vector_t, (K - 1)> &z_bar)
{
    const double dt = 1 / double(K - 1);
    using namespace boost::numeric::odeint;
    runge_kutta4<DiscretizationODE::state_type, double, DiscretizationODE::state_type, double, vector_space_algebra> stepper;

// #pragma omp parallel for schedule(dynamic)
    for (size_t k = 0; k < K - 1; k++)
    {
        DiscretizationODE::state_type V;
        V.setZero();
        V.col(0) = X.col(k);
        V.block<Model::state_dim_, Model::state_dim_>(0, 1).setIdentity();

        DiscretizationODE discretizationODE(U.col(k), U.col(k + 1), sigma, dt, model);
        integrate_adaptive(stepper, discretizationODE, V, 0., dt, dt / 10.);

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