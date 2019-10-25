#pragma once

#include "activeModel.hpp"

namespace scpp
{
op::SecondOrderConeProgram buildSCProblem(
    double &weight_time,
    double &weight_trust_region_time,
    double &weight_trust_region_trajectory,
    double &weight_virtual_control,
    TrajectoryData &td,
    DiscretizationData &dd);
}
