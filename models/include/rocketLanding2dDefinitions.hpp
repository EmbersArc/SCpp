#pragma once

#include <Eigen/Dense>

namespace rocket2d
{
    
enum rocket2d
{
    STATE_DIM_ = 6,
    INPUT_DIM_ = 2,
};

template <typename T>
void deg2rad(T &deg)
{
    deg *= M_PI / 180.;
}

} // namespace rocket2d