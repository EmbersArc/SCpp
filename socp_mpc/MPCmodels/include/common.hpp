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
    omegaMatrix << T(0.), -w(0), -w(1), -w(2),
        w(0), T(0.), w(2), -w(1),
        w(1), -w(2), T(0.), w(0),
        w(2), w(1), -w(0), T(0.);

    return omegaMatrix;
}

template <typename T>
Eigen::Matrix<T, 3, 3> omegaMatrixReduced(const Eigen::Matrix<T, 3, 1> &q)
{
    Eigen::Matrix<T, 3, 3> omegaMatrix;
    const T qw = sqrt(1. - q.squaredNorm());
    omegaMatrix << qw, -q(2), q(1),
        q(2), qw, -q(0),
        -q(1), q(0), qw;

    return omegaMatrix;
}

template <typename T>
Eigen::Matrix<T, 3, 3> EulerRotationMatrix(const Eigen::Matrix<T, 3, 1> &eta)
{
    const T phi = eta(0);
    const T theta = eta(1);
    const T psi = eta(2);

    Eigen::Matrix<T, 3, 3> R;
    R.row(0) << cos(theta) * cos(psi), -cos(theta) * sin(psi), sin(theta);
    R.row(1) << sin(phi) * sin(theta) * cos(psi) + cos(phi) * sin(psi),
        cos(phi) * cos(psi) - sin(phi) * sin(theta) * sin(psi), -sin(phi) * cos(theta);
    R.row(2) << -cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi),
        sin(phi) * cos(psi) + cos(phi) * sin(theta) * sin(psi), cos(phi) * cos(theta);

    return R;
}

template <typename T>
Eigen::Matrix<T, 3, 3> EulerRotationJacobian(const Eigen::Matrix<T, 3, 1> &eta)
{
    const T phi = eta(0);
    const T theta = eta(1);
    const T psi = eta(2);

    Eigen::Matrix<T, 3, 3> J;
    J.row(0) << cos(psi), -sin(psi), 0.;
    J.row(1) << cos(theta) * sin(psi), cos(theta) * cos(psi), 0;
    J.row(2) << -sin(theta) * cos(psi), sin(theta) * sin(psi), cos(theta);

    return 1. / cos(theta) * J;
}