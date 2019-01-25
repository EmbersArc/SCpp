#pragma once

#include "optimizationProblem.hpp"
#include "active_model.hpp"

namespace sc
{
op::SecondOrderConeProgram build_sc_SOCP(
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
    Model::state_vector_v_t &Sigma_bar,
    Model::state_vector_v_t &z_bar);
}

