#include "discretization.hpp"

namespace scpp::discretization
{

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

} // namespace scpp::discretization