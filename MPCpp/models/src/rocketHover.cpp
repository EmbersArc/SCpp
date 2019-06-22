#include "common.hpp"
#include "rocketHover.hpp"

using std::string;
using std::vector;

namespace rocketHover
{

RocketHover::RocketHover()
{
    p.loadFromFile();
}

void RocketHover::systemFlowMap(const state_vector_ad_t &x,
                                const input_vector_ad_t &u,
                                const param_vector_ad_t &par,
                                state_vector_ad_t &f)
{
    typedef scalar_ad_t T;

    // state variables
    auto v = x.segment<3>(3);
    auto q = x.segment<3>(6);
    auto w = x.segment<3>(9);

    auto R_I_B = Eigen::Quaternion<T>(sqrt(1. - q.squaredNorm()), q(0), q(1), q(2)).toRotationMatrix();
    auto J_B_inv = p.J_B.cast<T>().asDiagonal().inverse();
    auto g_I_ = p.g_I.cast<T>();
    auto r_T_B_ = p.r_T_B.cast<T>();

    f.segment<3>(0) << v;
    f.segment<3>(3) << 1. / T(p.m) * (R_I_B * u) + g_I_;
    f.segment<3>(6) << T(0.5) * omegaMatrixReduced<T>(q) * w;
    f.segment<3>(9) << J_B_inv * r_T_B_.cross(u) - w.cross(w);


    // // state variables
    // auto v = x.segment<3>(3);
    // auto eta = x.segment<3>(6);
    // auto w = x.segment<3>(9);

    // auto R_I_B = EulerRotationMatrix<T>(eta);
    // auto J_B_inv = p.J_B.cast<T>().asDiagonal().inverse();
    // auto g_I_ = p.g_I.cast<T>();
    // auto r_T_B_ = p.r_T_B.cast<T>();

    // f.segment<3>(0) << v;
    // f.segment<3>(3) << 1. / T(p.m) * (R_I_B * u) + g_I_;
    // f.segment<3>(6) << EulerRotationJacobian<T>(eta) * w;
    // f.segment<3>(9) << J_B_inv * r_T_B_.cross(u) - w.cross(w);
}

void RocketHover::getOperatingPoint(state_vector_t &x, input_vector_t &u)
{
    x << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    u << 0, 0, -p.g_I.z() * p.m;
}

void RocketHover::addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                            Eigen::MatrixXd &X0,
                                            Eigen::MatrixXd &U0)
{
    const size_t K = X0.cols();

    auto var = [&socp](const string &name, const vector<size_t> &indices = {}) { return socp.getVariable(name, indices); };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

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
        socp.addConstraint(op::norm2({(1.0) * var("X", {6, k}),
                                      (1.0) * var("X", {7, k})}) <= param_fn([this]() { return sqrt((1.0 - cos(p.theta_max)) / 2.); }));

        // Max Rotation Velocity
        socp.addConstraint(
            op::norm2({(1.0) * var("X", {9, k}),
                       (1.0) * var("X", {10, k}),
                       (1.0) * var("X", {11, k})}) <= param(p.w_B_max));
    }

    // Control Constraints
    for (size_t k = 0; k < K - 1; k++)
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
}

void RocketHover::nondimensionalize()
{
    // p.nondimensionalize();
}

void RocketHover::redimensionalize()
{
    // p.redimensionalize();
}

void RocketHover::Parameters::loadFromFile()
{
    ParameterServer param(fmt::format("../MPCpp/models/config/{}.info", getModelName()));

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

    x_init << r_init, v_init, eulerToQuaternion(rpy_init).tail<3>(), w_init;
    x_final << r_final, v_final, eulerToQuaternion(rpy_final).tail<3>(), w_final;
}

void RocketHover::Parameters::nondimensionalize()
{
    // m_scale = x_init(0);
    // r_scale = x_init.segment(0, 3).norm();

    // m /= m_scale;
    // r_T_B /= r_scale;
    // g_I /= r_scale;
    // J_B /= m_scale * r_scale * r_scale;

    // x_init.segment(0, 3) /= r_scale;
    // x_init.segment(3, 3) /= r_scale;

    // x_final.segment(0, 3) /= r_scale;
    // x_final.segment(3, 3) /= r_scale;

    // v_I_max /= r_scale;
    // T_min /= m_scale * r_scale;
    // T_max /= m_scale * r_scale;
}

void RocketHover::Parameters::redimensionalize()
{
    // m *= m_scale;
    // r_T_B *= r_scale;
    // g_I *= r_scale;
    // J_B *= m_scale * r_scale * r_scale;

    // x_init.segment(0, 3) *= r_scale;
    // x_init.segment(3, 3) *= r_scale;

    // x_final.segment(0, 3) *= r_scale;
    // x_final.segment(3, 3) *= r_scale;

    // v_I_max *= r_scale;
    // T_min *= m_scale * r_scale;
    // T_max *= m_scale * r_scale;
}

} // namespace rocketHover