#include "rocket2d.hpp"
#include "common.hpp"

namespace scpp::models
{

Rocket2d::Rocket2d() {}

void Rocket2d::systemFlowMap(const state_vector_ad_t &x,
                             const input_vector_ad_t &u,
                             const param_vector_ad_t &,
                             state_vector_ad_t &f)
{
    using T = scalar_ad_t;

    // state variables
    Eigen::Matrix<T, 2, 1> v(x(2), x(3));
    T eta = x(4);
    T w = x(5);

    // input variables
    T angle = u(0);
    T magnitude = u(1);
    Eigen::Matrix<T, 2, 1> T_B =
        Eigen::Rotation2D<T>(angle) * Eigen::Matrix<T, 2, 1>(0., magnitude);

    Eigen::Rotation2D<T> R_I_B(eta);

    f.segment<2>(0) << v;
    f.segment<2>(2) << 1. / T(p.m) * (R_I_B * T_B) + p.g_I.cast<T>();
    f(4) = w;
    f(5) = 1. / T(p.J_B) * (p.r_T_B.x() * T_B.y() - p.r_T_B.y() * T_B.x());
}

void Rocket2d::getOperatingPoint(state_vector_t &x, input_vector_t &u)
{
    x.setZero();
    u = -p.g_I * p.m;
}

void Rocket2d::addApplicationConstraints(op::SecondOrderConeProgram &socp,
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
        total_slack_variables = 4 * (X0.size() - 1); // four state constraints per timestep
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
        socp.addConstraint((1.0) * var("U", {0, X0.size() - 1}) == 0.0);
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
            socp.addConstraint(op::norm2({(1.0) * var("X", {2, k}),
                                          (1.0) * var("X", {3, k})}) <= rhs);
        }

        { // Max Tilt Angle
            op::AffineExpression rhs = param(p.theta_max);
            if (p.add_slack_variables)
            {
                rhs = rhs + 1.0 * var("epsilon", {slack_index++});
            }
            socp.addConstraint(op::norm2({(1.0) * var("X", {4, k})}) <= rhs);
        }

        { // Max Rotation Velocity
            op::AffineExpression rhs = param(p.w_B_max);
            if (p.add_slack_variables)
            {
                rhs = rhs + 1.0 * var("epsilon", {slack_index++});
            }
            socp.addConstraint(op::norm2({(1.0) * var("X", {5, k})}) <= rhs);
        }

        { // Glideslope
            op::AffineExpression rhs = param_fn([this]() { return std::tan(p.gamma_gs); }) * var("X", {1, k});
            if (p.add_slack_variables)
            {
                rhs = rhs + 1.0 * var("epsilon", {slack_index++});
            }
            socp.addConstraint(op::norm2({(1.0) * var("X", {0, k})}) <= rhs);
        }
    }
    if (p.add_slack_variables)
    {
        assert(slack_index == total_slack_variables);
    }

    // Control Constraints
    for (size_t k = 0; k < U0.size(); k++)
    {
        // Minimum Gimbal Angle
        socp.addConstraint((1.0) * var("U", {0, k}) + param(p.gimbal_max) >= (0.0));

        // Maximum Gimbal Angle
        socp.addConstraint((-1.0) * var("U", {0, k}) + param(p.gimbal_max) >= (0.0));

        // Minimum Thrust
        socp.addConstraint((1.0) * var("U", {1, k}) + -param(p.T_min) >= (0.0));

        // Maximum Thrust
        socp.addConstraint((-1.0) * var("U", {1, k}) + param(p.T_max) >= (0.0));
    }
}

void Rocket2d::nondimensionalize()
{
    p.nondimensionalize();
}

void Rocket2d::redimensionalize()
{
    p.redimensionalize();
}

void Rocket2d::nondimensionalizeTrajectory(trajectory_data_t &)
{
}

void Rocket2d::redimensionalizeTrajectory(trajectory_data_t &)
{
}

void Rocket2d::getInitializedTrajectory(trajectory_data_t &td)
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

void Rocket2d::loadParameters()
{
    p.loadFromFile(getParameterFolder() + "/model.info");
}

void Rocket2d::getNewModelParameters(param_vector_t &)
{
}

void Rocket2d::Parameters::loadFromFile(const std::string &path)
{
    ParameterServer param(path);

    param.loadMatrix("g_I", g_I);
    param.loadScalar("J_B", J_B);
    param.loadMatrix("r_T_B", r_T_B);

    Eigen::Vector2d r_init, v_init;
    Eigen::Vector2d r_final, v_final;
    double w_init, w_final;

    param.loadMatrix("r_init", r_init);
    param.loadMatrix("v_init", v_init);
    param.loadScalar("eta_init", eta_init);
    param.loadScalar("w_init", w_init);

    param.loadMatrix("r_final", r_final);
    param.loadMatrix("v_final", v_final);
    param.loadScalar("eta_final", eta_final);
    param.loadScalar("w_final", w_final);
    param.loadScalar("final_time", final_time);

    param.loadScalar("m", m);
    param.loadScalar("T_min", T_min);
    param.loadScalar("T_max", T_max);
    param.loadScalar("gamma_gs", gamma_gs);
    param.loadScalar("gimbal_max", gimbal_max);
    param.loadScalar("theta_max", theta_max);
    param.loadScalar("v_I_max", v_I_max);
    param.loadScalar("w_B_max", w_B_max);
    param.loadScalar("constrain_initial_final", constrain_initial_final);
    param.loadScalar("add_slack_variables", add_slack_variables);

    deg2rad(gimbal_max);
    deg2rad(theta_max);
    deg2rad(gamma_gs);
    deg2rad(w_B_max);
    deg2rad(w_init);
    deg2rad(w_final);
    deg2rad(eta_init);
    deg2rad(eta_final);

    x_init << r_init, v_init, eta_init, w_init;
    x_final << r_final, v_final, eta_final, w_final;
}

void Rocket2d::Parameters::nondimensionalize()
{
}

void Rocket2d::Parameters::redimensionalize()
{
}

} // namespace scpp::models