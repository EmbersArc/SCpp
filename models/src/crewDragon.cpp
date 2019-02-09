#include <random>

#include "common.hpp"
#include "crewDragon.hpp"
#include "parameterServer.hpp"

using std::string;
using std::vector;

namespace crewdragon
{

CrewDragon::CrewDragon()
{
    ParameterServer param(fmt::format("../models/config/{}.info", getModelName()));

    double I_sp;
    double m_init, m_dry;
    Eigen::Vector3d r_init, v_init, w_init;
    Eigen::Vector3d r_final, v_final;

    param.loadMatrix("g_I", g_I);
    param.loadMatrix("J_B", J_B);
    param.loadMatrix("r_T_B", r_T_B);
    param.loadScalar("m_init", m_init);
    param.loadMatrix("r_init", r_init);
    param.loadMatrix("v_init", v_init);
    param.loadMatrix("rpy_init", rpy_init);
    param.loadMatrix("w_init", w_init);
    param.loadScalar("m_dry", m_dry);
    param.loadMatrix("r_final", r_final);
    param.loadMatrix("v_final", v_final);
    param.loadScalar("final_time_guess", final_time_guess);
    param.loadScalar("T_min", T_min);
    param.loadScalar("T_max", T_max);
    param.loadScalar("I_sp", I_sp);
    param.loadScalar("theta_max", theta_max);
    param.loadScalar("gamma_gs", gamma_gs);
    param.loadScalar("w_B_max", w_B_max);

    deg2rad(gamma_gs);
    deg2rad(theta_max);
    deg2rad(w_init);
    deg2rad(w_B_max);
    deg2rad(rpy_init);

    alpha_m = 1. / (I_sp * fabs(g_I(2)));

    x_init << m_init, r_init, v_init, eulerToQuaternion(rpy_init), w_init;
    // randomizeInitialState();
    x_final << m_dry, r_final, v_final, 1., 0., 0., 0., 0, 0, 0;
}

void CrewDragon::systemFlowMap(
    const state_vector_ad_t &x,
    const input_vector_ad_t &u,
    state_vector_ad_t &f)
{
    typedef scalar_ad_t T;

    auto J_B_inv = J_B.cast<T>().asDiagonal().inverse();

    auto m = x(0);
    auto v_I = x.segment<3>(4);
    auto q_B_I = x.segment<4>(7);
    auto w_B = x.segment<3>(11);

    Eigen::AngleAxisd Rx(M_PI * 1. / 12., Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd Rz0(M_PI * 1. / 6., Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd Rz1(M_PI * 5. / 6., Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd Rz2(-M_PI * 5. / 6., Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd Rz3(-M_PI * 1. / 6., Eigen::Vector3d::UnitZ());

    Eigen::Matrix<T, 3, 1> u0 = (Rz0 * Rx * Eigen::Vector3d::UnitZ()).cast<T>() * u(0);
    Eigen::Matrix<T, 3, 1> u1 = (Rz1 * Rx * Eigen::Vector3d::UnitZ()).cast<T>() * u(1);
    Eigen::Matrix<T, 3, 1> u2 = (Rz2 * Rx * Eigen::Vector3d::UnitZ()).cast<T>() * u(2);
    Eigen::Matrix<T, 3, 1> u3 = (Rz3 * Rx * Eigen::Vector3d::UnitZ()).cast<T>() * u(3);

    f(0) = -T(alpha_m) * u.sum();
    f.segment(1, 3) << v_I;
    f.segment(4, 3) << 1. / m * dirCosineMatrix<T>(q_B_I).transpose() * (u0 + u1 + u2 + u3) + g_I.cast<T>();
    f.segment(7, 4) << T(0.5) * omegaMatrix<T>(w_B) * q_B_I;
    f.segment(11, 3) << J_B_inv * ((Rz0 * r_T_B).cast<T>().cross(u0) +
                                   (Rz1 * r_T_B).cast<T>().cross(u1) +
                                   (Rz2 * r_T_B).cast<T>().cross(u2) +
                                   (Rz3 * r_T_B).cast<T>().cross(u3)) -
                            w_B.cross(w_B);
}

void CrewDragon::initializeTrajectory(Eigen::MatrixXd &X,
                                      Eigen::MatrixXd &U)
{
    const size_t K = X.cols();

    for (size_t k = 0; k < K; k++)
    {
        const double alpha1 = double(K - k) / K;
        const double alpha2 = double(k) / K;

        // mass, position and linear velocity
        X(0, k) = alpha1 * x_init(0) + alpha2 * x_final(0);
        X.col(k).segment(1, 6) = alpha1 * x_init.segment(1, 6) + alpha2 * x_final.segment(1, 6);

        // do SLERP for quaternion
        Eigen::Quaterniond q0, q1;
        q0.w() = x_init(7);
        q0.vec() = x_init.segment(8, 3);
        q1.w() = x_final(7);
        q1.vec() << x_final.segment(8, 3);
        Eigen::Quaterniond qs = q0.slerp(alpha2, q1);
        X.col(k).segment(7, 4) << qs.w(), qs.vec();

        // angular velocity
        X.col(k).segment(11, 3) = alpha1 * x_init.segment(11, 3) + alpha2 * x_final.segment(11, 3);

        // input
        U.setConstant((T_max - T_min) / 2.);
    }
}

void CrewDragon::addApplicationConstraints(
    op::SecondOrderConeProgram &socp,
    Eigen::MatrixXd &X0,
    Eigen::MatrixXd &U0)
{
    const size_t K = X0.cols();

    auto var = [&socp](const string &name, const vector<size_t> &indices = {}) { return socp.getVariable(name, indices); };
    // auto param = [](double &param_value){ return op::Parameter(&param_value); };
    // auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    // Initial state
    for (size_t i = 0; i < STATE_DIM_; i++)
    {
        socp.addConstraint((-1.0) * var("X", {i, 0}) + (x_init(i)) == 0.0);
    }

    // Final State
    // mass and roll are free
    // socp.addConstraint((-1.0) * var("X", {0, K - 1}) + (x_final(0)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {1, K - 1}) + (x_final(1)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {2, K - 1}) + (x_final(2)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {3, K - 1}) + (x_final(3)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {4, K - 1}) + (x_final(4)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {5, K - 1}) + (x_final(5)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {6, K - 1}) + (x_final(6)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {7, K - 1}) + (x_final(7)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {8, K - 1}) + (x_final(8)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {9, K - 1}) + (x_final(9)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {10, K - 1}) + (x_final(10)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {11, K - 1}) + (x_final(11)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {12, K - 1}) + (x_final(12)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {13, K - 1}) + (x_final(13)) == 0.0);

    // State Constraints:
    for (size_t k = 0; k < K; k++)
    {
        // Mass
        //     x(0) >= m_dry
        //     for all k
        socp.addConstraint((1.0) * var("X", {0, k}) + (-x_final(0)) >= (0.0));

        // Max Tilt Angle
        // norm2([x(8), x(9)]) <= sqrt((1 - cos_theta_max) / 2)
        socp.addConstraint(op::norm2({(1.0) * var("X", {8, k}),
                                      (1.0) * var("X", {9, k})}) <= sqrt((1.0 - cos(theta_max)) / 2.));

        // Glide Slope
        socp.addConstraint(
            op::norm2({(1.0) * var("X", {1, k}),
                       (1.0) * var("X", {2, k})}) <= (1.0 / tan(gamma_gs)) * var("X", {3, k}));

        // Max Rotation Velocity
        socp.addConstraint(
            op::norm2({(1.0) * var("X", {11, k}),
                       (1.0) * var("X", {12, k}),
                       (1.0) * var("X", {13, k})}) <= (w_B_max));
    }

    // Control Constraints
    for (size_t k = 0; k < K; k++)
    {
        // Minimum and Maximum Thrust
        for (size_t i = 0; i < 4; i++)
        {
            socp.addConstraint((1.0) * var("U", {i, k}) + (-T_min) >= (0.0));
            socp.addConstraint((-1.0) * var("U", {i, k}) + (T_max) >= (0.0));
        }
    }
}

void CrewDragon::nondimensionalize()
{
    m_scale = x_init(0);
    r_scale = x_init.segment(1, 3).norm();

    alpha_m *= r_scale;
    r_T_B /= r_scale;
    g_I /= r_scale;
    J_B /= (m_scale * r_scale * r_scale);

    x_init(0) /= m_scale;
    x_init.segment(1, 3) /= r_scale;
    x_init.segment(4, 3) /= r_scale;

    x_final(0) /= m_scale;
    x_final.segment(1, 3) /= r_scale;
    x_final.segment(4, 3) /= r_scale;

    T_min /= m_scale * r_scale;
    T_max /= m_scale * r_scale;
}

void CrewDragon::redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                            Eigen::MatrixXd &U)
{
    const size_t K = X.cols();

    X.row(0) *= m_scale;
    X.block(1, 0, 6, K) *= r_scale;

    U *= m_scale * r_scale;
}

void CrewDragon::randomizeInitialState()
{
    std::mt19937 eng(time(0));
    auto dist = std::uniform_real_distribution<double>(-1., 1.);

    // mass
    x_init(0) *= 1.;

    // position
    x_init(1) *= dist(eng);
    x_init(2) *= dist(eng);
    x_init(3) *= 1.;

    // velocity
    x_init(4) *= dist(eng);
    x_init(5) *= dist(eng);
    x_init(6) *= 1. + 0.2 * dist(eng);

    // orientation
    double rx = dist(eng) * rpy_init.x();
    double ry = dist(eng) * rpy_init.y();
    double rz = rpy_init.z();
    Eigen::Vector3d euler(rx, ry, rz);
    x_init.segment(7, 4) << eulerToQuaternion(euler);
}

} // namespace crewdragon
