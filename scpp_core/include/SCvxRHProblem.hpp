#pragma once

#include "activeModel.hpp"

namespace scpp
{
op::SecondOrderConeProgram buildSCvxRHProblem(
    double &trust_region,
    double &weight_virtual_control,
    Model::state_vector_t &state_weights,
    Model::input_vector_t &input_weights,
    Model::state_vector_t &x_init,
    Model::state_vector_t &x_final,
    trajectory_data_t &td,
    discretization_data_t &dd);
}
