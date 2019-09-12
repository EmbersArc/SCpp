#include "rocketHover.hpp"
#include "common.hpp"

namespace rocketHover
{

RocketHover::RocketHover() {}

void RocketHover::systemFlowMap(const state_vector_ad_t &x,
                                const input_vector_ad_t &u,
                                const param_vector_ad_t &,
                                state_vector_ad_t &f)
{
    typedef scalar_ad_t T;

    // state variables
    Eigen::Matrix<T, 3, 1> v(x(3), x(4), x(5));
    Eigen::Matrix<T, 3, 1> eta(x(6), x(7), 0.);
    Eigen::Matrix<T, 3, 1> w(x(8), x(9), 0.);

    auto R_I_B = eulerToQuaternion(eta).toRotationMatrix();
    auto J_B_inv = p.J_B.cast<T>().asDiagonal().inverse();
    auto g_I_ = p.g_I.cast<T>();
    auto r_T_B_ = p.r_T_B.cast<T>();

    const Eigen::Matrix<T, 3, 1> thrust = u;

    f.segment<3>(0) << v;
    f.segment<3>(3) << 1. / T(p.m) * (R_I_B * thrust) + g_I_;
    f.segment<2>(6) << (rotationJacobian(eta) * w).head<2>();
    f.segment<2>(8) << (J_B_inv * (r_T_B_.cross(thrust)) - w.cross(w)).head<2>();
}

void RocketHover::getOperatingPoint(state_vector_t &x, input_vector_t &u)
{
    x << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    u << 0, 0, -p.g_I.z() * p.m;
}

void RocketHover::getStateWeights(state_vector_t &intermediate, state_vector_t &terminal)
{
    intermediate = p.state_weights_intermediate;
    terminal = p.state_weights_terminal;
}

void RocketHover::getInputWeights(input_vector_t &intermediate) { intermediate = p.input_weights; }

void RocketHover::addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                            Eigen::MatrixXd &X0,
                                            Eigen::MatrixXd &)
{
    const size_t K = X0.cols();

    auto var = [&socp](const std::string &name, const std::vector<size_t> &indices = {}) {
        return socp.getVariable(name, indices);
    };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    size_t total_slack_variables = 3 * (K - 1); // three state constraints per timestep
    socp.createTensorVariable("epsilon", {total_slack_variables});
    socp.createTensorVariable("epsilon_norm");
    std::vector<op::AffineExpression> norm2_terms;
    for (size_t i = 0; i < total_slack_variables; i++)
    {
        norm2_terms.push_back(1.0 * var("epsilon", {i}));
    }
    socp.addConstraint(op::norm2(norm2_terms) <= 1.0 * var("epsilon_norm"));
    socp.addMinimizationTerm(100. * var("epsilon_norm"));

    if (p.constrain_initial_final)
    {
        // Initial and final state
        for (size_t i = 0; i < STATE_DIM_; i++)
        {
            socp.addConstraint((-1.0) * var("X", {i, 0}) + param(p.x_init(i)) == 0.0);
            socp.addConstraint((-1.0) * var("X", {i, K - 1}) + param(p.x_final(i)) == 0.0);
        }
    }

    // State Constraints:
    size_t slack_index = 0;
    for (size_t k = 1; k < K; k++)
    {
        // Max Velocity
        socp.addConstraint(op::norm2({(1.0) * var("X", {3, k}),
                                      (1.0) * var("X", {4, k}),
                                      (1.0) * var("X", {5, k})}) <=
                           param(p.v_I_max) + 1.0 * var("epsilon", {slack_index++}));

        // Max Tilt Angle
        socp.addConstraint(op::norm2({(1.0) * var("X", {6, k}),
                                      (1.0) * var("X", {7, k})}) <=
                           param(p.theta_max) + 1.0 * var("epsilon", {slack_index++}));

        // Max Rotation Velocity
        socp.addConstraint(
            op::norm2({(1.0) * var("X", {8, k}),
                       (1.0) * var("X", {9, k})}) <=
            param(p.w_B_max) + 1.0 * var("epsilon", {slack_index++}));
    }
    assert(slack_index == total_slack_variables);

    // Control Constraints
    for (size_t k = 0; k < K - 1; k++)
    {
        // Simplified Minimum Thrust
        socp.addConstraint((1.0) * var("U", {2, k}) + param_fn([this]() { return -p.T_min; }) >= (0.0));

        // Maximum Thrust
        socp.addConstraint(op::norm2({(1.0) * var("U", {0, k}),
                                      (1.0) * var("U", {1, k}), (1.0) * var("U", {2, k})}) <=
                           param(p.T_max));

        // Maximum Gimbal Angle
        socp.addConstraint(op::norm2({(1.0) * var("U", {0, k}),
                                      (1.0) * var("U", {1, k})}) <=
                           param_fn([this]() { return std::tan(p.gimbal_max); }) * var("U", {2, k}));
    }
}

void RocketHover::nondimensionalize()
{
    // p.nondimensionalize();
}

void RocketHover::redimensionalize()
{
    // p.redimensionalize();
}

void RocketHover::getInitializedTrajectory(Eigen::MatrixXd &X,
                                           Eigen::MatrixXd &U)
{
    const size_t K = X.cols();

    for (size_t k = 0; k < K; k++)
    {
        const double alpha1 = double(K - k) / K;
        const double alpha2 = double(k) / K;

        // mass, position and linear velocity
        X.col(k) = alpha1 * p.x_init + alpha2 * p.x_final;

        // input
        U.setZero();
        U.row(2).setConstant(-p.g_I.z() * p.m);
    }
}

void RocketHover::loadParameters()
{
    p.loadFromFile(getParameterFolder() + "model.info");
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

    param.loadScalar("m", m);
    param.loadScalar("T_min", T_min);
    param.loadScalar("T_max", T_max);
    param.loadScalar("gimbal_max", gimbal_max);
    param.loadScalar("theta_max", theta_max);
    param.loadScalar("v_I_max", v_I_max);
    param.loadScalar("w_B_max", w_B_max);
    param.loadScalar("constrain_initial_final", constrain_initial_final);

    deg2rad(gimbal_max);
    deg2rad(theta_max);
    deg2rad(w_B_max);
    deg2rad(w_init);
    deg2rad(w_final);
    deg2rad(rpy_init);
    deg2rad(rpy_final);

    x_init << r_init, v_init, rpy_init.head<2>(), w_init.head<2>();
    x_final << r_final, v_final, rpy_final.head<2>(), w_final.head<2>();
}

void RocketHover::Parameters::nondimensionalize()
{
    m_scale = x_init(0);
    r_scale = x_init.segment(0, 3).norm();

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

} // namespace rocketHover