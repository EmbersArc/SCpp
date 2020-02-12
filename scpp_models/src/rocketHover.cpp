#include "rocketHover.hpp"
#include "common.hpp"

namespace scpp::models
{

void RocketHover::systemFlowMap(const state_vector_ad_t &x,
                                const input_vector_ad_t &u,
                                const param_vector_ad_t &par,
                                state_vector_ad_t &f)
{
    using T = scalar_ad_t;

    // state variables
    Eigen::Matrix<T, 3, 1> v(x(3), x(4), x(5));
    Eigen::Matrix<T, 3, 1> eta(x(6), x(7), T(0.));
    Eigen::Matrix<T, 3, 1> w(x(8), x(9), T(0.));

    auto R_I_B = eulerToQuaternionXYZ(eta).toRotationMatrix();

    auto m = par(0);
    auto g_I = par.segment<3>(1);
    auto J_B_inv = par.segment<3>(4).asDiagonal().inverse();
    auto r_T_B = par.segment<3>(7);
    // = 10 parameters

    f.segment<3>(0) << v;
    f.segment<3>(3) << 1. / m * (R_I_B * u) + g_I;
    f.segment<2>(6) << (rotationJacobianXYZ(eta) * w).head<2>();
    f.segment<2>(8) << (J_B_inv * (r_T_B.cross(u)) - w.cross(w)).head<2>();
}

void RocketHover::getOperatingPoint(state_vector_t &x, input_vector_t &u)
{
    x << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    u = -p.g_I * p.m;
}

void RocketHover::addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                            state_vector_v_t &,
                                            input_vector_v_t &)
{
    op::Variable v_epsilon;
    op::Variable v_X = socp.getVariable("X");
    op::Variable v_U = socp.getVariable("U");

    if (p.add_slack_variables)
    {
        v_epsilon = socp.createVariable("epsilon", 3, v_X.cols() - 1); // three state constraints
        op::Variable v_epsilon_norm = socp.createVariable("epsilon_norm");

        socp.addConstraint(op::norm2(v_epsilon) <= v_epsilon_norm);
        socp.addMinimizationTerm(op::Parameter(1000.) * v_epsilon_norm);
    }

    if (p.constrain_initial_final)
    {
        // Initial and final state
        socp.addConstraint(v_X.col(0) == op::Parameter(&p.x_init));
        socp.addConstraint(v_X.col(v_X.cols() - 1) == op::Parameter(&p.x_final));
        socp.addConstraint(v_U(0, v_U.cols() - 1) == op::Parameter(0.));
        socp.addConstraint(v_U(1, v_U.cols() - 1) == op::Parameter(0.));
    }

    // State Constraints:
    // Max Velocity
    op::SOCLhs lhs = op::norm2(v_X.block(3, 1, 3, v_X.cols() - 1), 0);
    if (p.add_slack_variables)
        lhs += v_epsilon.row(0);
    socp.addConstraint(lhs <= op::Parameter(&p.v_I_max));

    // Max Tilt Angle
    lhs = op::norm2(v_X.block(6, 1, 2, v_X.cols() - 1), 0);
    if (p.add_slack_variables)
        lhs += v_epsilon.row(1);
    socp.addConstraint(lhs <= op::Parameter(&p.theta_max));

    // Max Rotation Velocity
    lhs = op::norm2(v_X.block(8, 1, 2, v_X.cols() - 1), 0);
    if (p.add_slack_variables)
        lhs += v_epsilon.row(2);
    socp.addConstraint(lhs <= op::Parameter(&p.w_B_max));

    // Control Constraints
    // Simplified Minimum Thrust
    socp.addConstraint(v_U.row(2) >= op::Parameter(&p.T_min));

    // Maximum Thrust
    socp.addConstraint(op::norm2(v_U, 0) <= op::Parameter(&p.T_max));

    // Maximum Gimbal Angle
    socp.addConstraint(op::norm2(v_U.topRows(2), 0) <=
                       op::Parameter([this] { return std::tan(p.gimbal_max); }) * v_U.row(2));
}

void RocketHover::nondimensionalize()
{
    p.nondimensionalize();
}

void RocketHover::redimensionalize()
{
    p.redimensionalize();
}

void RocketHover::nondimensionalizeTrajectory(trajectory_data_t &td)
{
    for (auto &x : td.X)
    {
        x.segment<6>(0) /= p.r_scale;
    }
    for (auto &u : td.U)
    {
        u /= p.m_scale * p.r_scale;
    }
}

void RocketHover::redimensionalizeTrajectory(trajectory_data_t &td)
{
    for (auto &x : td.X)
    {
        x.segment<6>(0) *= p.r_scale;
    }
    for (auto &u : td.U)
    {
        u *= p.m_scale * p.r_scale;
    }
}

void RocketHover::getInitializedTrajectory(trajectory_data_t &td)
{
    for (size_t k = 0; k < td.n_X(); k++)
    {
        const double alpha1 = double(td.n_X() - k) / td.n_X();
        const double alpha2 = double(k) / td.n_X();

        td.X.at(k) = alpha1 * p.x_init + alpha2 * p.x_final;
    }

    std::fill(td.U.begin(), td.U.end(), -p.g_I * p.m);

    td.t = p.final_time;
}

void RocketHover::loadParameters()
{
    p.loadFromFile(getParameterFolder() + "/model.info");
}

void RocketHover::getNewModelParameters(param_vector_t &param)
{
    param << p.m, p.g_I, p.J_B, p.r_T_B;
}

void RocketHover::Parameters::loadFromFile(const std::string &path)
{
    ParameterServer param(path);

    param.loadMatrix("g_I", g_I);
    param.loadMatrix("J_B", J_B);
    param.loadMatrix("r_T_B", r_T_B);

    Eigen::Vector3d r_init, v_init, w_init;
    Eigen::Vector3d r_final, v_final, w_final;

    param.loadMatrix("r_init", r_init);
    param.loadMatrix("v_init", v_init);
    param.loadMatrix("rpy_init", rpy_init);
    param.loadMatrix("w_init", w_init);

    param.loadMatrix("r_final", r_final);
    param.loadMatrix("v_final", v_final);
    param.loadMatrix("rpy_final", rpy_final);
    param.loadMatrix("w_final", w_final);
    param.loadScalar("final_time", final_time);

    param.loadScalar("m", m);
    param.loadScalar("T_min", T_min);
    param.loadScalar("T_max", T_max);
    param.loadScalar("gimbal_max", gimbal_max);
    param.loadScalar("theta_max", theta_max);
    param.loadScalar("v_I_max", v_I_max);
    param.loadScalar("w_B_max", w_B_max);
    param.loadScalar("constrain_initial_final", constrain_initial_final);
    param.loadScalar("add_slack_variables", add_slack_variables);

    deg2rad(gimbal_max);
    deg2rad(theta_max);
    deg2rad(w_B_max);
    deg2rad(w_init);
    deg2rad(w_final);
    deg2rad(rpy_init);
    deg2rad(rpy_final);

    x_init << r_init, v_init, rpy_init.head<2>(), w_init.head<2>();
    x_final << r_final, v_final, rpy_final.head<2>(), w_final.head<2>();

    m_scale = x_init(0);
    r_scale = x_init.segment(0, 3).norm();
}

void RocketHover::Parameters::nondimensionalize()
{
    m /= m_scale;
    r_T_B /= r_scale;
    g_I /= r_scale;
    J_B /= m_scale * r_scale * r_scale;

    x_init.segment(0, 3) /= r_scale;
    x_init.segment(3, 3) /= r_scale;

    x_final.segment(0, 3) /= r_scale;
    x_final.segment(3, 3) /= r_scale;

    v_I_max /= r_scale;
    T_min /= m_scale * r_scale;
    T_max /= m_scale * r_scale;
}

void RocketHover::Parameters::redimensionalize()
{
    m *= m_scale;
    r_T_B *= r_scale;
    g_I *= r_scale;
    J_B *= m_scale * r_scale * r_scale;

    x_init.segment(0, 3) *= r_scale;
    x_init.segment(3, 3) *= r_scale;

    x_final.segment(0, 3) *= r_scale;
    x_final.segment(3, 3) *= r_scale;

    v_I_max *= r_scale;
    T_min *= m_scale * r_scale;
    T_max *= m_scale * r_scale;
}

} // namespace scpp::models