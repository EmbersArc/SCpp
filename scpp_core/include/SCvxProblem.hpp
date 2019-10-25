#pragma once

#include "activeModel.hpp"

namespace scpp
{

op::SecondOrderConeProgram buildSCvxProblem(
    double &trust_region,
    double &weight_virtual_control,
    Model::state_vector_v_t &X,
    Model::input_vector_v_t &U,
    DiscretizationData &dd);
}
