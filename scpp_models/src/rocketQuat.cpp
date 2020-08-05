#include "common.hpp"
#include "rocketQuat.hpp"

namespace scpp::models
{

void RocketQuat::systemFlowMap(const state_vector_ad_t &x,
                             const input_vector_ad_t &u,
                             const param_vector_ad_t &par,
                             state_vector_ad_t &f)
{
    using T = scalar_ad_t;

    auto alpha_m = par(0);
    auto g_I = par.segment<3>(1);
    auto J_B_inv = par.segment<3>(4).asDiagonal().inverse();
    auto r_T_B = par.segment<3>(7);
    // = 10 parameters

    // state variables
    auto m = x(0);
    auto v = x.segment<3>(4);
    auto q = x.segment<4>(7);
    auto w = x.segment<3>(11);

    auto thrust = u.head<3>();
    auto torque = Eigen::Matrix<T, 3, 1>(T(0.), T(0.), u(3));

    auto R_I_B = Eigen::Quaternion<T>(q(0), q(1), q(2), q(3))
                     .toRotationMatrix();

    f(0) = -alpha_m * thrust.norm();
    f.segment<3>(1) << v;
    f.segment<3>(4) << 1. / m * R_I_B * thrust + g_I;
    f.segment<4>(7) << T(0.5) * omegaMatrix<T>(w) * q;
    f.segment<3>(11) << J_B_inv * (r_T_B.cross(thrust) + torque) - w.cross(w);
}

void RocketQuat::getInitializedTrajectory(trajectory_data_t &td)
{
    for (size_t k = 0; k < td.n_X(); k++)
    {
        const double alpha1 = double(td.n_X() - k) / td.n_X();
        const double alpha2 = double(k) / td.n_X();

        // mass, position and linear velocity
        td.X.at(k)(0) = alpha1 * p.x_init(0) + alpha2 * p.x_final(0);
        td.X.at(k).segment(1, 6) = alpha1 * p.x_init.segment(1, 6) + alpha2 * p.x_final.segment(1, 6);

        // do SLERP for quaternion
        auto q0_ = p.x_init.segment<4>(7);
        auto q1_ = p.x_final.segment<4>(7);
        Eigen::Quaterniond q0(q0_(0), q0_(1), q0_(2), q0_(3));
        Eigen::Quaterniond q1(q1_(0), q1_(1), q1_(2), q1_(3));
        Eigen::Quaterniond qs = q0.slerp(alpha2, q1);
        td.X.at(k).segment<4>(7) << qs.w(), qs.vec();

        // angular velocity
        td.X.at(k).segment<3>(11) = alpha1 * p.x_init.segment<3>(11) + alpha2 * p.x_final.segment<3>(11);
    }

    for (auto &u : td.U)
    {
        u << 0., 0., (p.T_max - p.T_min) / 2., 0.;
    }

    td.t = p.final_time;
}

void RocketQuat::addApplicationConstraints(cvx::OptimizationProblem &socp,
                                           state_vector_v_t &,
                                           input_vector_v_t &U0)
{
    cvx::MatrixX v_X, v_U;
    socp.getVariable("X", v_X);
    socp.getVariable("U", v_U);

    // Initial state
    socp.addConstraint(cvx::equalTo(v_X.col(0), cvx::dynpar(p.x_init)));

    // Final State
    // mass and roll are free
    for (size_t i : {1, 2, 3,
                     4, 5, 6,
                     8, 9,
                     11, 12, 13})
    {
        socp.addConstraint(cvx::equalTo(v_X(i, v_X.cols() - 1), cvx::dynpar(p.x_final(i))));
    }

    // State Constraints:
    // Mass
    socp.addConstraint(cvx::greaterThan(v_X.row(0), cvx::dynpar(p.x_final(0))));

    // Glide Slope
    socp.addConstraint(cvx::lessThan(v_X.block(1, 0, 2, v_X.cols()).colwise().norm(),
                       cvx::dynpar(p_dyn.gs_const) * v_X.block(3, 0, 1, v_X.cols())));

    // Max Tilt Angle
    socp.addConstraint(cvx::lessThan(v_X.block(8, 0, 2, v_X.cols()).colwise().norm(),
                       cvx::dynpar(p_dyn.tilt_const)));

    // Max Rotation Velocity
    socp.addConstraint(cvx::lessThan(v_X.block(11, 0, 3, v_X.cols()).colwise().norm(),
                       cvx::dynpar(p.w_B_max)));

    // Control Constraints:
    // Final Input
    socp.addConstraint(cvx::equalTo(v_U.col(v_U.cols() - 1)(0), 0.));
    socp.addConstraint(cvx::equalTo(v_U.col(v_U.cols() - 1)(1), 0.));
    socp.addConstraint(cvx::equalTo(v_U.col(v_U.cols() - 1)(3), 0.));

    if (p.exact_minimum_thrust)
    {
        p_dyn.U0_ptr = &U0;
        p_dyn.thrust_const.resize(3, U0.size());

        // Linearized Minimum Thrust
        socp.addConstraint(cvx::greaterThan(cvx::dynpar(p_dyn.thrust_const).cwiseProduct(v_U.topRows(3)).colwise().sum(),
                           cvx::dynpar(p.T_min)));
    }
    else
    {
        // Simplified Minimum Thrust
        socp.addConstraint(cvx::greaterThan(v_U.row(2), cvx::dynpar(p.T_min)));
    }

    // Maximum Thrust
    socp.addConstraint(cvx::lessThan(v_U.topRows(3).colwise().norm(), cvx::dynpar(p.T_max)));

    // Maximum Gimbal Angle
    socp.addConstraint(cvx::lessThan(v_U.topRows(2).colwise().norm(),
                       cvx::dynpar(p_dyn.gimbal_const) * v_U.row(2)));

    if (p.enable_roll_control)
    {
        socp.addConstraint(cvx::box(-cvx::dynpar(p.t_max), v_U.row(3),cvx::dynpar(p.t_max)));
    }
    else
    {
        socp.addConstraint(cvx::equalTo(v_X.row(13) , 0.));
        socp.addConstraint(cvx::equalTo(v_U.row(3) , 0.));
    }
}

void RocketQuat::nondimensionalize()
{
    p.nondimensionalize();
}

void RocketQuat::redimensionalize()
{
    p.redimensionalize();
}

void RocketQuat::updateProblemParameters()
{
    p_dyn.gimbal_const = std::tan(p.gimbal_max);
    p_dyn.gs_const = std::tan(p.gamma_gs);
    p_dyn.tilt_const = std::sqrt((1. - std::cos(p.theta_max)) / 2.);

    for (size_t k = 0; k < size_t(p_dyn.thrust_const.cols()); k++)
    {
        p_dyn.thrust_const.col(k) = p_dyn.U0_ptr->at(k).head<3>().normalized();
    }
}

void RocketQuat::getNewModelParameters(param_vector_t &param)
{
    param << p.alpha_m, p.g_I, p.J_B, p.r_T_B;

    updateProblemParameters();
}

void RocketQuat::nondimensionalizeTrajectory(trajectory_data_t &td)
{
    for (auto &x : td.X)
    {
        x(0) /= p.m_scale;
        x.segment<6>(1) /= p.r_scale;
    }
    for (auto &u : td.U)
    {
        u.head<3>() /= p.m_scale * p.r_scale;
        u(3) /= p.m_scale * p.r_scale * p.r_scale;
    }
}

void RocketQuat::redimensionalizeTrajectory(trajectory_data_t &td)
{
    for (auto &x : td.X)
    {
        x(0) *= p.m_scale;
        x.segment<6>(1) *= p.r_scale;
    }
    for (auto &u : td.U)
    {
        u.head<3>() *= p.m_scale * p.r_scale;
        u(3) *= p.m_scale * p.r_scale * p.r_scale;
    }
}

void RocketQuat::Parameters::randomizeInitialState()
{
    // std::mt19937 eng(time(nullptr));
    // auto dist = std::uniform_real_distribution<double>(-1., 1.);

    // // mass
    // x_init(0) *= 1.;

    // // position
    // x_init(1) *= dist(eng);
    // x_init(2) *= dist(eng);
    // x_init(3) *= 1.;

    // // velocity
    // x_init(4) *= dist(eng);
    // x_init(5) *= dist(eng);
    // x_init(6) *= 1. + 0.2 * dist(eng);

    // // orientation
    // double rx = dist(eng) * rpy_init.x();
    // double ry = dist(eng) * rpy_init.y();
    // double rz = rpy_init.z();
    // Eigen::Vector3d euler(rx, ry, rz);
    // x_init.segment(7, 3) << quaternionToVector(eulerToQuaternionXYZ(euler));
}

void RocketQuat::loadParameters()
{
    p.loadFromFile(getParameterFolder() + "/model.info");
}

void RocketQuat::Parameters::loadFromFile(const std::string &path)
{
    ParameterServer param(path);

    bool random_initial_state;
    double I_sp;
    double m_init, m_dry;
    Eigen::Vector3d r_init, v_init, rpy_init, w_init;
    Eigen::Vector3d r_final, v_final, rpy_final, w_final;

    param.loadMatrix("g_I", g_I);
    param.loadMatrix("J_B", J_B);
    param.loadMatrix("r_T_B", r_T_B);
    param.loadScalar("m_init", m_init);
    param.loadMatrix("r_init", r_init);
    param.loadMatrix("v_init", v_init);
    param.loadMatrix("rpy_init", rpy_init);
    param.loadMatrix("w_init", w_init);
    param.loadMatrix("w_final", w_final);
    param.loadScalar("m_dry", m_dry);
    param.loadMatrix("r_final", r_final);
    param.loadMatrix("v_final", v_final);
    param.loadMatrix("rpy_final", rpy_final);
    param.loadScalar("T_min", T_min);
    param.loadScalar("T_max", T_max);
    param.loadScalar("t_max", t_max);
    param.loadScalar("I_sp", I_sp);
    param.loadScalar("gimbal_max", gimbal_max);
    param.loadScalar("theta_max", theta_max);
    param.loadScalar("gamma_gs", gamma_gs);
    param.loadScalar("w_B_max", w_B_max);
    param.loadScalar("random_initial_state", random_initial_state);
    param.loadScalar("final_time", final_time);
    param.loadScalar("exact_minimum_thrust", exact_minimum_thrust);
    param.loadScalar("enable_roll_control", enable_roll_control);

    deg2rad(gimbal_max);
    deg2rad(theta_max);
    deg2rad(gamma_gs);
    deg2rad(w_B_max);
    deg2rad(rpy_init);
    deg2rad(rpy_final);
    deg2rad(w_init);
    deg2rad(w_final);

    alpha_m = 1. / (I_sp * fabs(g_I(2)));

    const auto q_init = eulerToQuaternionXYZ(rpy_init);
    const auto q_final = eulerToQuaternionXYZ(rpy_final);
    x_init << m_init, r_init, v_init, q_init.w(), q_init.vec(), w_init;
    if (random_initial_state)
    {
        randomizeInitialState();
    }
    x_final << m_dry, r_final, v_final, q_final.w(), q_final.vec(), w_final;
}

void RocketQuat::Parameters::nondimensionalize()
{
    m_scale = x_init(0);
    r_scale = x_init.segment(1, 3).norm();

    alpha_m *= r_scale;
    r_T_B /= r_scale;
    g_I /= r_scale;
    J_B /= m_scale * r_scale * r_scale;

    x_init(0) /= m_scale;
    x_init.segment(1, 3) /= r_scale;
    x_init.segment(4, 3) /= r_scale;

    x_final(0) /= m_scale;
    x_final.segment(1, 3) /= r_scale;
    x_final.segment(4, 3) /= r_scale;

    T_min /= m_scale * r_scale;
    T_max /= m_scale * r_scale;
    t_max /= m_scale * r_scale * r_scale;
}

void RocketQuat::Parameters::redimensionalize()
{
    alpha_m /= r_scale;
    r_T_B *= r_scale;
    g_I *= r_scale;
    J_B *= m_scale * r_scale * r_scale;

    x_init(0) *= m_scale;
    x_init.segment(1, 3) *= r_scale;
    x_init.segment(4, 3) *= r_scale;

    x_final(0) *= m_scale;
    x_final.segment(1, 3) *= r_scale;
    x_final.segment(4, 3) *= r_scale;

    T_min *= m_scale * r_scale;
    T_max *= m_scale * r_scale;
    t_max *= m_scale * r_scale * r_scale;
}

} // namespace scpp::models