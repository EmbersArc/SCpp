#include "simulation.hpp"

const Eigen::Matrix<double, 3, 2> I3x2 = Eigen::Matrix<double, 3, 2>::Identity();
const Eigen::Matrix<double, 2, 3> I2x3 = Eigen::Matrix<double, 2, 3>::Identity();

Eigen::Vector3d toEulerAngle(const Eigen::Quaterniond &q)
{
    const Eigen::Matrix3d R = q.toRotationMatrix();
    const double phi = atan2(-R(1, 2), R(2, 2));
    const double theta = asin(R(0, 2));
    const double psi = -atan2(-R(0, 1), R(0, 0));

    return Eigen::Vector3d(phi, theta, psi);
}

void coneProject(Eigen::Vector3d &v, double angle, double cutoff, double height)
{
    v.z() = std::clamp(v.z(), cutoff, height);

    const double current_angle = fabs(atan(v.head<2>().norm() / v.z()));
    if (current_angle > angle)
    {
        // scale x and y components to be within gimbal constraint
        v.head<2>() *= tan(angle) * v.z() / v.head<2>().norm();
    }
}

inline Eigen::Vector3d vee(const Eigen::Matrix3d &S) { return Eigen::Vector3d(S(2, 1), S(0, 2), S(1, 0)); }

Eigen::Vector3d customPositionController(std::shared_ptr<Model> model, Model::state_vector_t &x_init,
                                         Model::state_vector_t &x_final)
{
    const Eigen::Vector3d position_gain_(1, 1, 1);
    const Eigen::Vector3d velocity_gain_(1, 1, 1);
    const Eigen::Vector3d attitude_gain_(20, 20, 1);
    const Eigen::Vector3d angular_rate_gain_(1, 1, 1);
    const double kMaxVel_ = 2.;
    const double yaw_des = 0.;
    // const double kMaxAngle_ = 0.2;
    // const double kMinThrust_ = 5;
    // const double kMaxThrust_ = 20;

    const Eigen::Vector3d r = x_init.segment<3>(0);
    const Eigen::Vector3d v = x_init.segment<3>(3);
    const Eigen::Vector3d q_red = x_init.segment<3>(6);
    const Eigen::Vector3d omega = x_init.segment<3>(9);

    const Eigen::Vector3d lambda_ref = x_final.segment<3>(0);
    // const Eigen::Vector3d sigma_ref = x_final.segment<3>(3);
    const Eigen::Vector3d a_ref(0, 0, 0);
    const Eigen::Vector3d q_ref_red = x_final.segment<3>(6);
    const Eigen::Vector3d omega_ref = x_final.segment<3>(9);

    // reduced quaternion to full quaternion
    Eigen::Quaterniond q;
    q.w() = std::sqrt(1. - q_red.dot(q_red));
    q.vec() << q_red;
    Eigen::Quaterniond q_ref;
    q_ref.w() = std::sqrt(1. - q_ref_red.dot(q_ref_red));
    q_ref.vec() << q_ref_red;

    const Eigen::Matrix3d R = q.toRotationMatrix();
    const Eigen::Vector3d E_z = Eigen::Vector3d::UnitZ();

    const double km_ = model->p.m;
    const Eigen::Vector3d kI_ = model->p.J_B;
    const double kg_ = -model->p.g_I.z();
    const double kl_ = -model->p.r_T_B.z();
    const double kc_ = kI_.x() / (km_ * kl_);

    const Eigen::Vector3d lambda = r + kc_ * (R * E_z - E_z);
    const Eigen::Vector3d sigma = v + omega.cross(kc_ * R * E_z);

    // desired position
    const Eigen::Vector3d e_x = lambda_ref - lambda;

    // desired velocity
    Eigen::Vector3d sigma_des = position_gain_.asDiagonal() * e_x;
    if (sigma_des.norm() > kMaxVel_)
    {
        sigma_des = sigma_des.normalized() * kMaxVel_;
    }

    // velocity error
    const Eigen::Vector3d e_v = sigma_des - sigma;

    // desired force
    const Eigen::Vector3d f_des = velocity_gain_.asDiagonal() * e_v + km_ * kg_ * E_z + km_ * a_ref;

    // desired orientation
    const double yaw = yaw_des;
    Eigen::Vector3d b_d(cos(yaw), sin(yaw), 0);

    Eigen::Vector3d b_c_3;
    b_c_3 << f_des.x(), f_des.y(), f_des.z();
    b_c_3.normalize();
    const Eigen::Vector3d b_c_2 = b_c_3.cross(b_d);
    Eigen::Matrix3d R_c;
    R_c.col(0) = b_c_2.cross(b_c_3);
    R_c.col(1) = b_c_2;
    R_c.col(2) = b_c_3;

    const Eigen::Vector3d angle_error = vee(0.5 * (R.transpose() * R_c - R_c.transpose() * R));
    const Eigen::Vector3d angular_velocity_error = omega_ref - omega;

    // desired angular velocity
    const Eigen::Vector3d omega_des = attitude_gain_.asDiagonal() * angle_error +
                                      angular_rate_gain_.asDiagonal() * angular_velocity_error +
                                      omega.cross(kI_.asDiagonal() * omega);

    const Eigen::Vector3d t_des = kI_.asDiagonal().inverse() * omega_des;

    // desired lateral force
    const Eigen::Vector2d TyTx = I2x3 * t_des / kl_;

    // desired longitudinal force
    const double Tz = f_des.dot(R * E_z);

    Eigen::Vector3d control_vector(-TyTx[1], TyTx[0], Tz);

    control_vector.z() = std::max(0., control_vector.z());
    double rx = std::atan(-control_vector.y() / control_vector.z());
    double ry = std::atan(control_vector.x() / control_vector.z());
    rx = std::clamp(rx, -model->p.gimbal_max, model->p.gimbal_max);
    ry = std::clamp(ry, -model->p.gimbal_max, model->p.gimbal_max);
    const double thrust_magnitude = std::max(0., std::min(control_vector.head<3>().norm(), model->p.T_max));
    control_vector = Eigen::AngleAxisd(rx, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()) *
                     (Eigen::Vector3d::UnitZ() * thrust_magnitude);

    return control_vector;
}

Eigen::Vector3d woodPositionController(std::shared_ptr<Model> model, Model::state_vector_t &x_init,
                                       Model::state_vector_t &x_final)
{
    const Eigen::Vector3d lambda_des = x_final.segment<3>(0);

    const Eigen::Vector3d r = x_init.segment<3>(0);
    const Eigen::Vector3d v = x_init.segment<3>(3);
    const Eigen::Vector3d q_red = x_init.segment<3>(6);
    const Eigen::Vector3d omega = x_init.segment<3>(9);

    Eigen::Quaterniond q;
    q.w() = std::sqrt(1. - q_red.dot(q_red));
    q.vec() << q_red;
    const Eigen::Matrix3d R = q.toRotationMatrix();
    const Eigen::Vector3d E_z = Eigen::Vector3d::UnitZ();

    const double kI1_ = model->p.J_B.x();
    const double km_ = model->p.m;
    const double kg_ = -model->p.g_I.z();
    const double kl_ = -model->p.r_T_B.z();
    const double kc_ = kI1_ / (km_ * kl_);
    const Eigen::Vector3d lambda = (r + kc_ * (R * E_z - E_z)) - lambda_des;
    const Eigen::Vector3d sigma = v + omega.cross(kc_ * R * E_z);
    const Eigen::Vector3d eta = toEulerAngle(q);

    const double k_1_ = 1;
    const double k_2_ = 3;
    const double k_3_ = 1;
    const double k_4_ = 3;

    //! begin generated
    const double lambda0 = lambda(0, 0);
    const double lambda1 = lambda(1, 0);
    const double lambda2 = lambda(2, 0);
    const double sigma0 = sigma(0, 0);
    const double sigma1 = sigma(1, 0);
    const double sigma2 = sigma(2, 0);
    const double eta0 = eta(0, 0);
    const double eta1 = eta(1, 0);
    const double eta2 = eta(2, 0);
    const double omega0 = omega(0, 0);
    const double omega1 = omega(1, 0);
    const double omega2 = omega(2, 0);
    const double x0 = sin(eta2);
    const double x1 = cos(eta1);
    const double x2 = k_3_ * k_4_ + 1;
    const double x3 = k_1_ * k_2_ + 1;
    const double x4 = lambda1 * x3;
    const double x5 = -k_1_ - k_2_;
    const double x6 = sigma1 * x5;
    const double x7 = sigma2 * x5;
    const double x8 = kg_ - lambda2 * x3 + x7;
    const double x9 = cos(eta2);
    const double x10 = 1.0 * omega0;
    const double x11 = x0 * x10;
    const double x12 = 1.0 * omega1;
    const double x13 = x12 * x9;
    const double x14 = x11 + x13;
    const double x15 = sin(eta1);
    const double x16 = 1.0 * x15;
    const double x17 = x14 * x16 / pow(x1, 2);
    const double x18 = 1.0 / x1;
    const double x19 = x10 * x18 * x9;
    const double x20 = x0 * x12 * x18;
    const double x21 = 1.0 * omega2 - x15 * x19 + x15 * x20;
    const double x22 = x0 * x21;
    const double x23 = 1.0 * x18;
    const double x24 = x21 * x9;
    const double x25 = k_3_ + k_4_;
    const double x26 = -x4 + x6;
    const double x27 = pow(x26, 2) + pow(x8, 2);
    const double x28 = sigma0 * x5;
    const double x29 = -lambda0 * x3 + x28;
    const double x30 = x27 + pow(x29, 2);
    const double x31 = sqrt(x30);
    const double x32 = 1.0 * x31;
    const double x33 = x3 * x32;
    const double x34 = sin(eta0);
    const double x35 = x1 * x34;
    const double x36 = x33 * x35;
    const double x37 = x36 + x6;
    const double x38 = x37 * x8;
    const double x39 = cos(eta0);
    const double x40 = x1 * x39;
    const double x41 = -x3 * (-kg_ + x32 * x40);
    const double x42 = x41 + x7;
    const double x43 = x26 * x42;
    const double x44 = x26 * x37 + x42 * x8;
    const double x45 = x0 * x39;
    const double x46 = x34 * x9;
    const double x47 = x12 * x31;
    const double x48 = x39 * x9;
    const double x49 = x0 * x34;
    const double x50 = x10 * x31;
    const double x51 = -x15 * x33;
    const double x52 = x28 + x51;
    const double x53 = x29 * x52 + x44;
    const double x54 = x53 / x31;
    const double x55 = 1.0 * x54;
    const double x56 = x36 + x5 * (-x35 * x55 + x47 * (x15 * x46 + x45) - x50 * (-x15 * x49 + x48));
    const double x57 = x41 + x5 * (x40 * x55 + x47 * (-x15 * x48 + x49) - x50 * (x15 * x45 + x46));
    const double x58 = 1.0 * kI1_ / kl_;
    const double x59 =
        x1 * x58 *
        (-omega0 * (x17 * x9 - x22 * x23) - omega1 * (-x0 * x17 - x23 * x24) - x2 * (eta0 - atan2(x4 - x6, x8)) -
         x25 * (x19 - x20 - (-x38 + x43) / x27) + (x27 * (x26 * x57 - x56 * x8) + 2 * x44 * (x38 - x43)) / pow(x27, 2));
    const double x60 = sqrt(x27);
    const double x61 = x27 * x52;
    const double x62 = x29 * x44;
    const double x63 = -x61 + x62;
    const double x64 = x1 * x31;
    const double x65 =
        x58 * (-x10 * x24 + x12 * x22 - x2 * (eta1 - atan2(x29, x60)) - x25 * (x14 - (x61 - x62) / (x30 * x60)) +
               (x27 * x30 *
                    (x27 * (x5 * (x11 * x64 + x13 * x64 + x16 * x54) + x51) -
                     x29 * (x26 * x56 + pow(x37, 2) + pow(x42, 2) + x57 * x8) + x44 * x52) +
                2 * x27 * x53 * x63 + x30 * x44 * x63) /
                   (pow(x27, 3.0 / 2.0) * pow(x30, 2)));

    Eigen::Vector3d control_vector;
    control_vector.setZero();

    control_vector(0, 0) = x0 * x59 - x65 * x9;
    control_vector(1, 0) = x0 * x65 + x59 * x9;
    control_vector(2, 0) = km_ * x31;

    //! end generated

    control_vector.z() = std::max(0.01, control_vector.z());

    double rx = std::atan(-control_vector.y() / control_vector.z());
    double ry = std::atan(control_vector.x() / control_vector.z());

    rx = std::clamp(rx, -model->p.gimbal_max, model->p.gimbal_max);
    ry = std::clamp(ry, -model->p.gimbal_max, model->p.gimbal_max);

    const double thrust_magnitude = std::max(0., std::min(control_vector.head<3>().norm(), model->p.T_max));
    control_vector = Eigen::AngleAxisd(rx, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()) *
                     (Eigen::Vector3d::UnitZ() * thrust_magnitude);

    return control_vector;
}
