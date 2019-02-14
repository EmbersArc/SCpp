#pragma once

#include <Eigen/Dense>

template <typename T>
void deg2rad(T &deg)
{
    deg *= M_PI / 180.;
}

template <typename T>
Eigen::Matrix<T, 4, 1> eulerToQuaternion(const Eigen::Matrix<T, 3, 1> &rpy)
{
    Eigen::Quaternion<T> q = Eigen::AngleAxis(rpy.x(), Eigen::Matrix<T, 3, 1>::UnitX()) *
                             Eigen::AngleAxis(rpy.y(), Eigen::Matrix<T, 3, 1>::UnitY()) *
                             Eigen::AngleAxis(rpy.z(), Eigen::Matrix<T, 3, 1>::UnitZ());
    Eigen::Matrix<T, 4, 1> quat;
    quat << q.w(), q.vec();

    return quat;
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