#pragma once

#include "optimizationProblem.hpp"
#include "activeModel.hpp"

namespace mpc
{
op::SecondOrderConeProgram buildSCOP(
    Model::ptr_t model,
    Model::state_vector_v_t &X,
    Model::input_vector_v_t &U,
    Model::state_vector_t &x_init,
    Model::state_vector_t &x_final,
    Model::state_vector_t &state_weights_intermediate,
    Model::state_vector_t &state_weights_terminal,
    Model::input_vector_t &input_weights,
    Model::state_matrix_t &A,
    Model::control_matrix_t &B,
    Model::state_vector_t &z,
    bool constant_dynamics = false,
    bool intermediate_cost_active = true);
}
