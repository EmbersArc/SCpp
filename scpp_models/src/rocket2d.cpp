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
                                         state_vector_v_t &,
                                         input_vector_v_t &)
{
    op::Variable v_X = socp.getVariable("X");
    op::Variable v_U = socp.getVariable("U");

    if (p.constrain_initial_final)
    {
        // Initial and final state
        socp.addConstraint(-v_X.col(0) + op::Parameter(&p.x_init) == 0.);
        socp.addConstraint(-v_X.col(v_X.cols() - 1) + op::Parameter(&p.x_final) == 0.);
        socp.addConstraint(op::Parameter(1.0) * v_U(0, v_U.cols() - 1) == 0.);
    }

    // State Constraints:
    // Glideslope
    socp.addConstraint(op::Norm2(v_X.row(0), 0) <=
                       op::Parameter([this]() { return std::tan(p.gamma_gs); }) * v_X.row(1));
    // Max Velocity
    socp.addConstraint(op::Norm2(v_X.block(2, 0, 2, v_X.cols()), 0) <=
                       op::Parameter(&p.v_I_max));
    // Max Tilt Angle
    socp.addConstraint(op::Norm2(v_X.row(4), 0) <=
                       op::Parameter(&p.theta_max));
    // Max Rotation Velocity
    socp.addConstraint(op::Norm2(v_X.row(5), 0) <=
                       op::Parameter(&p.w_B_max));

    // Control Constraints
    auto some_vector = Eigen::MatrixXd(1, v_U.cols());
    some_vector.setOnes();
    // Minimum Gimbal Angle
    socp.addConstraint(v_U.row(0) + op::Parameter(&p.gimbal_max) * op::Parameter(some_vector) >= (0.0));
    // Maximum Gimbal Angle
    socp.addConstraint(-v_U.row(0) + op::Parameter(&p.gimbal_max) * op::Parameter(some_vector) >= (0.0));
    // Minimum Thrust
    socp.addConstraint(v_U.row(1) + op::Parameter(&p.T_min) * op::Parameter(-some_vector) >= (0.0));
    // Maximum Thrust
    socp.addConstraint(-v_U.row(1) + op::Parameter(&p.T_max) * op::Parameter(some_vector) >= (0.0));
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