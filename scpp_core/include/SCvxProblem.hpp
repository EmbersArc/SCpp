#pragma once

#include "activeModel.hpp"

namespace scpp
{

op::SecondOrderConeProgram buildSCvxProblem(
    double &trust_region,
    double &weight_virtual_control,
    TrajectoryData &td,
    DiscretizationData &dd);
}
