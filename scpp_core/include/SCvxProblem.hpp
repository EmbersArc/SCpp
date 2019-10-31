#pragma once

#include "activeModel.hpp"

namespace scpp
{

op::SecondOrderConeProgram buildSCvxProblem(
    double &trust_region,
    double &weight_virtual_control,
    trajectory_data_t &td,
    discretization_data_t &dd);
}
