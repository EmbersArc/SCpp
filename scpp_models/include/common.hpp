#pragma once

#include <Eigen/Dense>

namespace scpp::models
{

template <typename T>
void deg2rad(T &deg)
{
    deg *= M_PI / 180.;
}

template <typename T>
void rad2deg(T &rad)
{
    rad *= 180. / M_PI;
}

template <typename T>
Eigen::Matrix<T, 4, 1> quaternionToVector(const Eigen::Quaternion<T> &q)
{
    Eigen::Matrix<T, 4, 1> q_vec;
    q_vec << q.w(), q.vec();
    return q_vec;
}

// sequence x-y'-z'' is XYZ = Z''Y'X
template <typename T>
Eigen::Quaternion<T> eulerToQuaternionXYZ(const Eigen::Matrix<T, 3, 1> &eta)
{
    Eigen::Quaternion<T> q;
    q = Eigen::AngleAxis<T>(eta.x(), Eigen::Matrix<T, 3, 1>::UnitX()) *
        Eigen::AngleAxis<T>(eta.y(), Eigen::Matrix<T, 3, 1>::UnitY()) *
        Eigen::AngleAxis<T>(eta.z(), Eigen::Matrix<T, 3, 1>::UnitZ());

    return q;
}

template <typename T>
Eigen::Matrix<T, 3, 3> eulerRotationMatrixXY(const Eigen::Matrix<T, 2, 1> &eta)
{
    const T phi = eta.x();
    const T theta = eta.y();

    Eigen::Matrix<T, 3, 3> M;
    M.row(0) << cos(theta), T(0.), sin(theta);
    M.row(1) << sin(theta) * sin(phi), cos(phi), -sin(phi) * cos(theta);
    M.row(2) << -sin(theta) * cos(phi), sin(phi), cos(phi) * cos(theta);

    return M;
}

// sequence x-y-z is ZYX
template <typename T>
Eigen::Quaternion<T> eulerToQuaternionZYX(const Eigen::Matrix<T, 3, 1> &eta)
{
    Eigen::Quaternion<T> q;
    q = Eigen::AngleAxis<T>(eta.z(), Eigen::Matrix<T, 3, 1>::UnitZ()) *
        Eigen::AngleAxis<T>(eta.y(), Eigen::Matrix<T, 3, 1>::UnitY()) *
        Eigen::AngleAxis<T>(eta.x(), Eigen::Matrix<T, 3, 1>::UnitX());
    return q;
}

template <typename T>
Eigen::Matrix<T, 3, 1> quaternionToEulerXYZ(const Eigen::Quaternion<T> &q)
{
    const Eigen::Matrix<T, 3, 3> R = q.toRotationMatrix();
    const T phi = atan2(-R(1, 2), R(2, 2));
    const T theta = asin(R(0, 2));
    const T psi = atan2(-R(0, 1), R(0, 0));
    return Eigen::Matrix<T, 3, 1>(phi, theta, psi);
}

template <typename T>
Eigen::Matrix<T, 3, 1> quaternionToEulerZYX(const Eigen::Quaternion<T> &q)
{
    const Eigen::Matrix<T, 3, 3> R = q.toRotationMatrix();
    const T phi = atan2(R(1, 0), R(0, 0));
    const T theta = asin(-R(2, 0));
    const T psi = atan2(R(2, 1), R(2, 2));
    return Eigen::Matrix<T, 3, 1>(psi, theta, phi);
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
Eigen::Matrix<T, 3, 3> rotationJacobianXYZ(const Eigen::Matrix<T, 3, 1> &eta)
{
    // const T phi = eta.x();
    const T theta = eta.y();
    const T psi = eta.z();

    Eigen::Matrix<T, 3, 3> M;
    M.row(0) << cos(psi), -sin(psi), T(0.);
    M.row(1) << cos(theta) * sin(psi), cos(theta) * cos(psi), T(0.);
    M.row(2) << -sin(theta) * cos(psi), sin(theta) * sin(psi), cos(theta);

    return M / cos(theta);
}

template <typename T>
Eigen::Matrix<T, 2, 2> rotationJacobianXY(const Eigen::Matrix<T, 2, 1> &eta)
{
    const T theta = eta.y();

    Eigen::Matrix<T, 2, 2> M;
    M.row(0) << cos(theta), sin(theta);
    M.row(1) << -sin(theta), cos(theta);

    return M;
}

template <typename T>
Eigen::Matrix<T, 4, 4> omegaMatrix(const Eigen::Matrix<T, 3, 1> &w)
{
    Eigen::Matrix<T, 4, 4> omega;
    omega << T(0.), -w(0), -w(1), -w(2),
        w(0), T(0.), w(2), -w(1),
        w(1), -w(2), T(0.), w(0),
        w(2), w(1), -w(0), T(0.);

    return omega;
}

template <typename T>
Eigen::Matrix<T, 3, 3> omegaMatrixReduced(const Eigen::Matrix<T, 3, 1> &q)
{
    Eigen::Matrix<T, 3, 3> omega;
    const T qw = sqrt(1. - q.squaredNorm());
    omega << qw, -q(2), q(1),
        q(2), qw, -q(0),
        -q(1), q(0), qw;

    return omega;
}

} // namespace scpp::models