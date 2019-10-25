#pragma once

#include "activeModel.hpp"

namespace scpp
{
op::SecondOrderConeProgram buildSCProblem(
    double &weight_time,
    double &weight_trust_region_time,
    double &weight_trust_region_trajectory,
    double &weight_virtual_control,
    Model::state_vector_v_t &X,
    Model::input_vector_v_t &U,
    double &sigma,
    DiscretizationData &dd);
}
