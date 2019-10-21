#pragma once

#include "optimizationProblem.hpp"
#include "activeModel.hpp"

namespace sc
{
op::SecondOrderConeProgram buildSCOP(
    Model::ptr_t model,
    double &trust_region,
    double &weight_virtual_control,
    Model::state_vector_v_t &X,
    Model::input_vector_v_t &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &z_bar);
}
