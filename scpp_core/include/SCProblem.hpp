#pragma once

#include "activeModel.hpp"

namespace scpp
{
cvx::OptimizationProblem buildSCProblem(
    double &weight_time,
    double &weight_trust_region_time,
    double &weight_trust_region_trajectory,
    double &weight_virtual_control,
    trajectory_data_t &td,
    discretization_data_t &dd);
}
