#include "common.hpp"
#include "rocket3d.hpp"

using std::string;
using std::vector;

namespace rocket3d
{

Rocket3D::Rocket3D()
{
    p.loadFromFile();
}

void Rocket3D::systemFlowMap(const state_vector_ad_t &x,
                             const input_vector_ad_t &u,
                             const param_vector_ad_t &p,
                             state_vector_ad_t &f)
{
    typedef scalar_ad_t T;

    auto m = p(0);
    auto J_B_inv = p.segment<3>(1).asDiagonal().inverse();
    auto g_I_ = p.segment<3>(4);
    auto r_T_B_ = p.segment<3>(7);
    // = 10 parameters

    // state variables
    auto v = x.segment<3>(3);
    auto q = x.segment<4>(6);
    auto w = x.segment<3>(10);

    auto R_I_B = Eigen::Quaternion<T>(q(0), q(1), q(2), q(3)).toRotationMatrix();

    f.segment<3>(0) << v;
    f.segment<3>(3) << 1. / m * R_I_B * u + g_I_;
    f.segment<4>(6) << T(0.5) * omegaMatrix<T>(w) * q;
    f.segment<3>(10) << J_B_inv * r_T_B_.cross(u) - w.cross(w);
}

void Rocket3D::getInitializedTrajectory(Eigen::MatrixXd &X,
                                        Eigen::MatrixXd &U)
{
    const size_t K = X.cols();

    for (size_t k = 0; k < K; k++)
    {
        const double a = double(K - k) / K;
        const double b = double(k) / K;

        // position and linear velocity
        X.col(k).segment(0, 3) = a * p.x_init.segment(0, 3) + b * p.x_final.segment(0, 3);
        X.col(k).segment(3, 3) = (p.x_final.segment(0, 3) - p.x_init.segment(0, 3)).normalized() * p.v_I_max / 2;

        // do SLERP for quaternion
        Eigen::Quaterniond q0(p.x_init(6), p.x_init(7), p.x_init(8), p.x_init(9));
        Eigen::Quaterniond q1(p.x_final(6), p.x_final(7), p.x_final(8), p.x_final(9));
        Eigen::Quaterniond qs = q0.slerp(b, q1);
        X.col(k).segment(6, 4) << qs.w(), qs.vec();

        // angular velocity
        X.col(k).segment(10, 3) = a * p.x_init.segment(10, 3) + b * p.x_final.segment(10, 3);

        // input
        U.setZero();
        U.row(2).setConstant((p.T_max - p.T_min) / 2.);
    }
}

void Rocket3D::addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                         Eigen::MatrixXd &X0,
                                         Eigen::MatrixXd &U0)
{
    const size_t K = X0.cols();

    auto var = [&socp](const string &name, const vector<size_t> &indices = {}) { return socp.getVariable(name, indices); };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    // Initial state
    for (size_t i = 0; i < state_dim; i++)
    {
        socp.addConstraint((-1.0) * var("X", {i, 0}) + param(p.x_init(i)) == 0.0);
    }
    // // Final state
    // socp.addConstraint((-1.0) * var("X", {0, K - 1}) + param(p.x_final(0)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {1, K - 1}) + param(p.x_final(1)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {2, K - 1}) + param(p.x_final(2)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {3, K - 1}) + param(p.x_final(3)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {4, K - 1}) + param(p.x_final(4)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {5, K - 1}) + param(p.x_final(5)) == 0.0);
    // // socp.addConstraint((-1.0) * var("X", {6, K - 1}) + param(p.x_final(6)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {7, K - 1}) + param(p.x_final(7)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {8, K - 1}) + param(p.x_final(8)) == 0.0);
    // // socp.addConstraint((-1.0) * var("X", {9, K - 1}) + param(p.x_final(9)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {10, K - 1}) + param(p.x_final(10)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {11, K - 1}) + param(p.x_final(11)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {12, K - 1}) + param(p.x_final(12)) == 0.0);

    // // Final state
    // socp.createTensorVariable("error", {state_dim}); // error minimization term
    // socp.createTensorVariable("error_norm");         // error minimization term
    // std::vector<op::AffineExpression> norm_args;
    // for (size_t i = 0; i < state_dim; i++)
    // {
    //     socp.addConstraint((-1.0) * var("X", {i, K - 1}) + param(p.x_final(i)) + 1.0 * var("error", {i}) == 0.0);
    //     norm_args.push_back(param(p.state_weights(i)) * var("error", {i}));
    // }
    // socp.addConstraint(op::norm2(norm_args) <= 1.0 * var("error_norm"));
    // socp.addMinimizationTerm(1.0 * var("error_norm"));

    // State Constraints:
    for (size_t k = 0; k < K; k++)
    {
        // Max Velocity
        socp.addConstraint(
            op::norm2({(1.0) * var("X", {3, k}),
                       (1.0) * var("X", {4, k}),
                       (1.0) * var("X", {5, k})}) <= param(p.v_I_max));

        // Max Tilt Angle
        // norm2([x(7), x(8)]) <= sqrt((1 - cos_theta_max) / 2)
        socp.addConstraint(op::norm2({(1.0) * var("X", {7, k}),
                                      (1.0) * var("X", {8, k})}) <= param_fn([this]() { return sqrt((1.0 - cos(p.theta_max)) / 2.); }));

        // Max Rotation Velocity
        socp.addConstraint(
            op::norm2({(1.0) * var("X", {10, k}),
                       (1.0) * var("X", {11, k}),
                       (1.0) * var("X", {12, k})}) <= param(p.w_B_max));
    }

    // Control Constraints
    for (size_t k = 0; k < K; k++)
    {
        // Simplified Minimum Thrust
        socp.addConstraint((1.0) * var("U", {2, k}) + param_fn([this]() { return -p.T_min; }) >= (0.0));

        // Maximum Thrust
        socp.addConstraint(
            op::norm2({(1.0) * var("U", {0, k}),
                       (1.0) * var("U", {1, k}),
                       (1.0) * var("U", {2, k})}) <= param(p.T_max));

        // Maximum Gimbal Angle
        socp.addConstraint(
            op::norm2({(1.0) * var("U", {0, k}),
                       (1.0) * var("U", {1, k})}) <= param_fn([this]() { return tan(p.gimbal_max); }) * var("U", {2, k}));
    }

    /**
     * Build error cost
     *
     */
    std::vector<op::AffineExpression> error_norm2_args;
    for (size_t k = 0; k < K; k++)
    {
        for (size_t i = 0; i < state_dim; i++)
        {
            op::AffineExpression ex = param_fn([this, i, k]() { return -1.0 * p.state_weights(i) * p.x_final(i); }) + param(p.state_weights(i)) * var("X", {i, k});
            error_norm2_args.push_back(ex);
        }
    }
    socp.createTensorVariable("error_cost"); // error minimization term
    socp.addConstraint(op::norm2(error_norm2_args) <= (1.0) * var("error_cost"));
    socp.addMinimizationTerm(1.0 * var("error_cost"));
}

void Rocket3D::nondimensionalize()
{
    p.nondimensionalize();
}

void Rocket3D::redimensionalize()
{
    p.redimensionalize();
}

void Rocket3D::nondimensionalizeTrajectory(Eigen::MatrixXd &X,
                                           Eigen::MatrixXd &U)
{
    p.nondimensionalizeTrajectory(X, U);
}

void Rocket3D::redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                          Eigen::MatrixXd &U)
{
    p.redimensionalizeTrajectory(X, U);
}

void Rocket3D::getNewModelParameters(param_vector_t &param)
{
    param << p.m, p.J_B, p.g_I, p.r_T_B;
}

void Rocket3D::Parameters::loadFromFile()
{
    ParameterServer param(fmt::format("../SCpp/models/config/{}.info", getModelName()));

    Eigen::Vector3d r_init, v_init, rpy_init, w_init;
    Eigen::Vector3d r_final, v_final, rpy_final, w_final;

    param.loadMatrix("g_I", g_I);
    param.loadMatrix("J_B", J_B);
    param.loadMatrix("r_T_B", r_T_B);

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

    deg2rad(gimbal_max);
    deg2rad(theta_max);
    deg2rad(w_B_max);
    deg2rad(w_init);
    deg2rad(w_final);
    deg2rad(rpy_init);
    deg2rad(rpy_final);

    param.loadMatrix("state_weights", state_weights);

    x_init << r_init, v_init, eulerToQuaternion(rpy_init), w_init;
    x_final << r_final, v_final, eulerToQuaternion(rpy_final), w_final;
}

void Rocket3D::Parameters::nondimensionalize()
{
    m_scale = m;
    r_scale = x_init.segment(0, 3).norm();

    m /= m_scale;
    J_B /= m_scale * r_scale * r_scale;
    g_I /= r_scale;
    r_T_B /= r_scale;

    x_init.segment(0, 3) /= r_scale;
    x_init.segment(3, 3) /= r_scale;

    x_final.segment(0, 3) /= r_scale;
    x_final.segment(3, 3) /= r_scale;

    v_I_max /= r_scale;
    T_min /= m_scale * r_scale;
    T_max /= m_scale * r_scale;
}

void Rocket3D::Parameters::redimensionalize()
{
    m *= m_scale;
    J_B *= m_scale * r_scale * r_scale;
    g_I *= r_scale;
    r_T_B *= r_scale;

    x_init.segment(0, 3) *= r_scale;
    x_init.segment(3, 3) *= r_scale;

    x_final.segment(0, 3) *= r_scale;
    x_final.segment(3, 3) *= r_scale;

    v_I_max *= r_scale;
    T_min *= m_scale * r_scale;
    T_max *= m_scale * r_scale;
}

void Rocket3D::Parameters::nondimensionalizeTrajectory(Eigen::MatrixXd &X, Eigen::MatrixXd &U) const
{
    const size_t K = X.cols();

    X.block(0, 0, 6, K) /= r_scale;

    U /= m_scale * r_scale;
}

void Rocket3D::Parameters::redimensionalizeTrajectory(Eigen::MatrixXd &X, Eigen::MatrixXd &U) const
{
    const size_t K = X.cols();

    X.block(0, 0, 6, K) *= r_scale;

    U *= m_scale * r_scale;
}

} // namespace rocket3d