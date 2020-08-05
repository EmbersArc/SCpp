#pragma once

#include "activeModel.hpp"

namespace scpp
{

cvx::OptimizationProblem buildSCvxProblem(
    double &trust_region,
    double &weight_virtual_control,
    trajectory_data_t &td,
    discretization_data_t &dd);
}
