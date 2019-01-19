#pragma once

#include <Eigen/Dense>

namespace rocket3d
{
enum rocket3d
{
    STATE_DIM_ = 14,
    INPUT_DIM_ = 3,
};

template <typename T>
void deg2rad(T &deg)
{
    deg *= M_PI / 180.;
}

template <typename T>
Eigen::Matrix<T, 3, 3> skew(const Eigen::Matrix<T, 3, 1> &v)
{
    Eigen::Matrix<T, 3, 3> skewMatrix;
    skewMatrix << T(0), -v(2), v(1),
        v(2), T(0), -v(0),
        -v(1), v(0), T(0);
    return skewMatrix;
}

template <typename T>
Eigen::Matrix<T, 3, 3> dirCosineMatrix(const Eigen::Matrix<T, 4, 1> &q)
{
    Eigen::Matrix<T, 3, 3> dirCosineMatrix;
    dirCosineMatrix << 1 - 2 * (q(2) * q(2) + q(3) * q(3)), 2 * (q(1) * q(2) + q(0) * q(3)), 2 * (q(1) * q(3) - q(0) * q(2)),
        2 * (q(1) * q(2) - q(0) * q(3)), 1 - 2 * (q(1) * q(1) + q(3) * q(3)), 2 * (q(2) * q(3) + q(0) * q(1)),
        2 * (q(1) * q(3) + q(0) * q(2)), 2 * (q(2) * q(3) - q(0) * q(1)), 1 - 2 * (q(1) * q(1) + q(2) * q(2));

    return dirCosineMatrix;
}

template <typename T>
Eigen::Matrix<T, 4, 4> omegaMatrix(const Eigen::Matrix<T, 3, 1> &w)
{
    Eigen::Matrix<T, 4, 4> omegaMatrix;
    omegaMatrix << T(0), -w(0), -w(1), -w(2),
        w(0), T(0), w(2), -w(1),
        w(1), -w(2), T(0), w(0),
        w(2), w(1), -w(0), T(0);

    return omegaMatrix;
}

} // namespace rocket3d