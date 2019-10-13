#include "common.hpp"
#include "starship.hpp"

using std::string;
using std::vector;

namespace scpp::models
{

Starship::Starship() {}

void Starship::systemFlowMap(const state_vector_ad_t &x,
                             const input_vector_ad_t &u,
                             const param_vector_ad_t &par,
                             state_vector_ad_t &f)
{
    typedef scalar_ad_t T;

    auto alpha_m_ = par(0);
    auto g_I_ = par.segment<3>(1);
    auto J_B_inv = par.segment<3>(4).asDiagonal().inverse();
    auto r_T_B_ = par.segment<3>(7);
    // = 10 parameters

    // state variables
    auto m = x(0);
    auto v = x.segment<3>(4);
    auto q = x.segment<4>(7);
    auto w = x.segment<3>(11);

    auto R_I_B = Eigen::Quaternion<T>(q(0), q(1), q(2), q(3)).toRotationMatrix();

    f(0) = -alpha_m_ * u.norm();
    f.segment(1, 3) << v;
    f.segment(4, 3) << 1. / m * R_I_B * u + g_I_;
    f.segment(7, 4) << T(0.5) * omegaMatrix<T>(w) * q;
    f.segment(11, 3) << J_B_inv * r_T_B_.cross(u) - w.cross(w);
}

void Starship::getInitializedTrajectory(state_vector_v_t &X,
                                        input_vector_v_t &U,
                                        double &t)
{
    for (size_t k = 0; k < X.size(); k++)
    {
        const double alpha1 = double(X.size() - k) / X.size();
        const double alpha2 = double(k) / X.size();

        // mass, position and linear velocity
        X.at(k)(0) = alpha1 * p.x_init(0) + alpha2 * p.x_final(0);
        X.at(k).segment(1, 6) = alpha1 * p.x_init.segment(1, 6) + alpha2 * p.x_final.segment(1, 6);

        // do SLERP for quaternion
        Eigen::Quaterniond q0(p.x_init(7), p.x_init(8), p.x_init(9), p.x_init(10));
        Eigen::Quaterniond q1(p.x_final(7), p.x_final(8), p.x_final(9), p.x_final(10));
        Eigen::Quaterniond qs = q0.slerp(alpha2, q1);
        X.at(k).segment(7, 4) << qs.w(), qs.vec();

        // angular velocity
        X.at(k).segment(11, 3) = alpha1 * p.x_init.segment(11, 3) + alpha2 * p.x_final.segment(11, 3);
    }

    for (size_t k = 0; k < U.size(); k++)
    {
        // input
        U.at(k) << 0, 0, (p.T_max - p.T_min) / 2.;
    }
    t = p.final_time;
}

void Starship::addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                         state_vector_v_t &X0,
                                         input_vector_v_t &U0)
{
    const size_t K = X0.size();

    auto var = [&socp](const string &name, const vector<size_t> &indices = {}) { return socp.getVariable(name, indices); };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    // Initial state
    for (size_t i = 0; i < STATE_DIM_; i++)
    {
        socp.addConstraint((-1.0) * var("X", {i, 0}) + param(p.x_init(i)) == 0.0);
    }

    // Final State
    // mass and roll are free
    for (size_t i : {1, 2, 3, 4, 5, 6, 8, 9, 11, 12})
    {
        socp.addConstraint((-1.0) * var("X", {i, K - 1}) + param(p.x_final(i)) == 0.0);
    }

    // Final Input
    socp.addConstraint((1.0) * var("U", {0, K - 1}) == (0.0));
    socp.addConstraint((1.0) * var("U", {1, K - 1}) == (0.0));

    // State Constraints:
    for (size_t k = 0; k < K; k++)
    {
        // Mass
        //     x(0) >= m_dry
        //     for all k
        socp.addConstraint(
            (1.0) * var("X", {0, k}) + param_fn([this]() { return -p.x_final(0); }) >= (0.0));

        // Glide Slope
        socp.addConstraint(
            op::norm2({(1.0) * var("X", {1, k}),
                       (1.0) * var("X", {2, k})}) <= param_fn([this]() { return tan(p.gamma_gs); }) * var("X", {3, k}));

        // Max Tilt Angle
        // norm2([x(8), x(9)]) <= sqrt((1 - cos_theta_max) / 2)
        socp.addConstraint(
            op::norm2({(1.0) * var("X", {8, k}),
                       (1.0) * var("X", {9, k})}) <=
            param_fn([this]() { return sqrt((1.0 - cos(p.theta_max)) / 2.); }));

        // Max Rotation Velocity
        socp.addConstraint(
            op::norm2({(1.0) * var("X", {11, k}),
                       (1.0) * var("X", {12, k}),
                       (1.0) * var("X", {13, k})}) <= param(p.w_B_max));
    }

    // Control Constraints
    for (size_t k = 0; k < K; k++)
    {
        // Linearized Minimum Thrust
        op::AffineExpression lhs;
        for (size_t i = 0; i < INPUT_DIM_; i++)
        {
            lhs = lhs + param_fn([&U0, i, k]() { return (U0.at(k)(i) / sqrt(U0.at(k)(0) * U0.at(k)(0) + U0.at(k)(1) * U0.at(k)(1) + U0.at(k)(2) * U0.at(k)(2))); }) * var("U", {i, k});
        }
        socp.addConstraint(lhs + param_fn([this]() { return -p.T_min; }) >= (0.0));

        // Simplified Minimum Thrust
        // socp.addConstraint((1.0) * var("U", {2, k}) + param_fn([this]() { return -p.T_min; }) >= (0.0));

        // Maximum Thrust
        socp.addConstraint(
            op::norm2({(1.0) * var("U", {0, k}),
                       (1.0) * var("U", {1, k}),
                       (1.0) * var("U", {2, k})}) <= param(p.T_max));

        // Maximum Gimbal Angle
        socp.addConstraint(
            op::norm2({(1.0) * var("U", {0, k}),
                       (1.0) * var("U", {1, k})}) <=
            param_fn([this]() { return tan(p.gimbal_max); }) * var("U", {2, k}));
    }
}

void Starship::nondimensionalize()
{
    p.nondimensionalize();
}

void Starship::redimensionalize()
{
    p.redimensionalize();
}

void Starship::getNewModelParameters(param_vector_t &param)
{
    param << p.alpha_m, p.g_I, p.J_B, p.r_T_B;
}

void Starship::nondimensionalizeTrajectory(state_vector_v_t &X,
                                           input_vector_v_t &U)
{
    for (size_t k = 0; k < X.size(); k++)
    {
        X[k](0) /= p.m_scale;
        X[k].segment<6>(1) /= p.r_scale;
    }
    for (size_t k = 0; k < U.size(); k++)
    {
        U[k] /= p.m_scale * p.r_scale;
    }
}

void Starship::redimensionalizeTrajectory(state_vector_v_t &X,
                                          input_vector_v_t &U)
{
    for (size_t k = 0; k < X.size(); k++)
    {
        X[k](0) *= p.m_scale;
        X[k].segment<6>(1) *= p.r_scale;
    }
    for (size_t k = 0; k < U.size(); k++)
    {
        U[k] *= p.m_scale * p.r_scale;
    }
}

void Starship::Parameters::randomizeInitialState()
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
    x_init.segment(7, 4) << quaternionToVector(eulerToQuaternion(euler));
}

void Starship::loadParameters()
{
    p.loadFromFile(getParameterFolder() + "/model.info");
}

void Starship::Parameters::loadFromFile(const std::string &path)
{
    ParameterServer param(path);

    bool random_initial_state;
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
    param.loadScalar("T_min", T_min);
    param.loadScalar("T_max", T_max);
    param.loadScalar("I_sp", I_sp);
    param.loadScalar("gimbal_max", gimbal_max);
    param.loadScalar("theta_max", theta_max);
    param.loadScalar("gamma_gs", gamma_gs);
    param.loadScalar("w_B_max", w_B_max);
    param.loadScalar("random_initial_state", random_initial_state);
    param.loadScalar("final_time", final_time);

    deg2rad(gimbal_max);
    deg2rad(theta_max);
    deg2rad(gamma_gs);
    deg2rad(w_init);
    deg2rad(w_B_max);
    deg2rad(rpy_init);

    alpha_m = 1. / (I_sp * fabs(g_I(2)));

    Eigen::Vector4d q_init = quaternionToVector(eulerToQuaternion(rpy_init));
    x_init << m_init, r_init, v_init, q_init, w_init;
    if (random_initial_state)
    {
        randomizeInitialState();
    }
    x_final << m_dry, r_final, v_final, 1., 0., 0., 0., 0, 0, 0;
}

void Starship::Parameters::nondimensionalize()
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
}

void Starship::Parameters::redimensionalize()
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
}

} // namespace scpp::models