#pragma once

#include "optimizationProblem.hpp"
#include "activeModel.hpp"

namespace sc
{
op::SecondOrderConeProgram buildSCOP(
    Model &model,
    double &weight_trust_region_time,
    double &weight_trust_region_trajectory,
    double &weight_virtual_control,
    Eigen::MatrixXd &X,
    Eigen::MatrixXd &U,
    double &sigma,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &S_bar,
    Model::state_vector_v_t &z_bar);
}

