#pragma once

#include <Eigen/Dense>
#include <iostream>

template <typename T>
void deg2rad(T &deg)
{
    deg *= M_PI / 180.;
}

template <typename T>
Eigen::Quaternion<T> eulerToQuaternion(const Eigen::Matrix<T, 3, 1> &eta)
{
    Eigen::Quaternion<T> q;
    q = Eigen::AngleAxis<T>(eta.x(), Eigen::Matrix<T, 3, 1>::UnitX()) *
        Eigen::AngleAxis<T>(eta.y(), Eigen::Matrix<T, 3, 1>::UnitY()) *
        Eigen::AngleAxis<T>(eta.z(), Eigen::Matrix<T, 3, 1>::UnitZ());
    return q;
}

template <typename T>
Eigen::Matrix<T, 4, 1> quaternionToVector(const Eigen::Quaternion<T> &q)
{
    Eigen::Matrix<T, 4, 1> q_vec;
    q_vec << q.w(), q.vec();
    return q_vec;
}

template <typename T>
Eigen::Quaternion<T> vectorToQuaternion(const Eigen::Matrix<T, 3, 1> &v)
{
    return Eigen::Quaternion<T>(std::sqrt(1. - v.squaredNorm()), v.x(), v.y(), v.z());
}

template <typename T>
Eigen::Quaternion<T> vectorToQuaternion(const Eigen::Matrix<T, 4, 1> &v)
{
    return Eigen::Quaternion<T>(v(3), v(0), v(1), v(2));
}

template <typename T>
Eigen::Matrix<T, 3, 1> quaternionToEulerRPY(const Eigen::Quaternion<T> &q)
{
    const Eigen::Matrix<T, 3, 3> R = q.toRotationMatrix();
    const T phi = atan2(-R(1, 2), R(2, 2));
    const T theta = asin(R(0, 2));
    const T psi = atan2(-R(0, 1), R(0, 0));

    return Eigen::Matrix<T, 3, 1>(phi, theta, psi);
}

template <typename T>
Eigen::Matrix<T, 3, 1> quaternionToEulerYPR(const Eigen::Quaternion<T> &q)
{
    const Eigen::Matrix<T, 3, 3> R = q.toRotationMatrix();
    const T phi = atan2(R(2, 1), R(2, 2));
    const T theta = asin(-R(2, 0));
    const T psi = atan2(R(1, 0), R(0, 0));
    return Eigen::Matrix<T, 3, 1>(phi, theta, psi);
}

template <typename T>
Eigen::Matrix<T, 3, 3> rotationJacobian(const Eigen::Matrix<T, 3, 1> &eta)
{
    const T phi = eta.x();
    const T theta = eta.y();
    const T psi = eta.z();

    Eigen::Matrix<T, 3, 3> M;
    M.row(0) << cos(psi), -sin(psi), 0;
    M.row(1) << cos(theta) * sin(psi), cos(theta) * cos(psi), 0;
    M.row(2) << -sin(theta) * cos(psi), sin(theta) * sin(psi), cos(theta);

    return M / cos(theta);
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
