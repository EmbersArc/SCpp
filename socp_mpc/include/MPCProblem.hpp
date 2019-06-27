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
    const Model::state_matrix_t &A,
    const Model::control_matrix_t &B,
    const Model::state_vector_t &z);
}
