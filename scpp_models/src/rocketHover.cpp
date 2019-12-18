#include "rocketHover.hpp"
#include "common.hpp"

namespace scpp::models
{

RocketHover::RocketHover() {}

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
                                            state_vector_v_t &X0,
                                            input_vector_v_t &U0)
{
    auto var = [&socp](const std::string &name, const std::vector<size_t> &indices = {}) {
        return socp.getVariable(name, indices);
    };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    size_t total_slack_variables;
    if (p.add_slack_variables)
    {
        total_slack_variables = 3 * (X0.size() - 1); // three state constraints per timestep
        socp.createTensorVariable("epsilon", {total_slack_variables});
        socp.createTensorVariable("epsilon_norm");

        std::vector<op::AffineExpression> norm2_terms;
        for (size_t i = 0; i < total_slack_variables; i++)
        {
            norm2_terms.push_back(1.0 * var("epsilon", {i}));
        }
        socp.addConstraint(op::norm2(norm2_terms) <= 1.0 * var("epsilon_norm"));
        socp.addMinimizationTerm(1000. * var("epsilon_norm"));
    }

    if (p.constrain_initial_final)
    {
        // Initial and final state
        for (size_t i = 0; i < STATE_DIM_; i++)
        {
            socp.addConstraint((-1.0) * var("X", {i, 0}) + param(p.x_init(i)) == 0.0);
            socp.addConstraint((-1.0) * var("X", {i, X0.size() - 1}) + param(p.x_final(i)) == 0.0);
        }
    }

    // State Constraints:
    size_t slack_index = 0;
    for (size_t k = 1; k < X0.size(); k++)
    {
        { // Max Velocity
            op::AffineExpression rhs = param(p.v_I_max);
            if (p.add_slack_variables)
            {
                rhs = rhs + 1.0 * var("epsilon", {slack_index++});
            }
            socp.addConstraint(op::norm2({(1.0) * var("X", {3, k}),
                                          (1.0) * var("X", {4, k}),
                                          (1.0) * var("X", {5, k})}) <= rhs);
        }

        { // Max Tilt Angle
            op::AffineExpression rhs = param(p.theta_max);
            if (p.add_slack_variables)
            {
                rhs = rhs + 1.0 * var("epsilon", {slack_index++});
            }
            socp.addConstraint(op::norm2({(1.0) * var("X", {6, k}),
                                          (1.0) * var("X", {7, k})}) <= rhs);
        }

        { // Max Rotation Velocity
            op::AffineExpression rhs = param(p.w_B_max);
            if (p.add_slack_variables)
            {
                rhs = rhs + 1.0 * var("epsilon", {slack_index++});
            }
            socp.addConstraint(
                op::norm2({(1.0) * var("X", {8, k}),
                           (1.0) * var("X", {9, k})}) <= rhs);
        }
    }
    if (p.add_slack_variables)
    {
        assert(slack_index == total_slack_variables);
    }

    // Control Constraints
    for (size_t k = 0; k < U0.size(); k++)
    {
        // Simplified Minimum Thrust
        socp.addConstraint((1.0) * var("U", {2, k}) + -param(p.T_min) >= (0.0));

        // Maximum Thrust
        socp.addConstraint(op::norm2({(1.0) * var("U", {0, k}),
                                      (1.0) * var("U", {1, k}),
                                      (1.0) * var("U", {2, k})}) <=
                           param(p.T_max));

        // Maximum Gimbal Angle
        socp.addConstraint(op::norm2({(1.0) * var("U", {0, k}),
                                      (1.0) * var("U", {1, k})}) <=
                           param_fn([this]() { return std::tan(p.gimbal_max); }) * var("U", {2, k}));
    }
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