#pragma once

#include "optimizationProblem.hpp"
#include "activeModel.hpp"

namespace mpc
{
op::SecondOrderConeProgram buildSCOP(
    Model &model,
    Eigen::MatrixXd &X,
    Eigen::MatrixXd &U,
    Model::state_vector_t &x_init,
    Model::state_vector_t &x_final,
    Model::state_matrix_t &A,
    Model::control_matrix_t &B,
    Model::state_vector_t &z);
}
