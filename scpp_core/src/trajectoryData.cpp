#include "trajectoryData.hpp"

void TrajectoryData::initialize(size_t K, bool interpolate_input, bool free_final_time)
{
    A.resize(K - 1);
    B.resize(K - 1);
    if (interpolate_input)
    {
        C.resize(K - 1);
    }
    if (free_final_time)
    {
        s.resize(K - 1);
    }
    z.resize(K - 1);
}

bool TrajectoryData::interpolatedInput() const
{
    return not C.empty();
}

bool TrajectoryData::variableTime() const
{
    return not s.empty();
}
