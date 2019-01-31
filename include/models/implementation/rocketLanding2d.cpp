#include "models/rocketLanding2d.hpp"
#include "parameterServer.hpp"

namespace rocket2d
{

RocketLanding2D::RocketLanding2D()
{
    ParameterServer param(format("../include/models/config/{}.info", getModelName()));

    param.loadScalar("m", m);
    param.loadScalar("g", g);
    param.loadScalar("r_T", r_T);
    param.loadScalar("I", I);
    param.loadScalar("T_min", T_min);
    param.loadScalar("T_max", T_max);
    param.loadScalar("gimbal_max", gimbal_max);
    param.loadScalar("theta_max", theta_max);
    param.loadScalar("gamma_gs", gamma_gs);

    param.loadMatrix("x_init", x_init);
    param.loadMatrix("x_final", x_final);

    deg2rad(x_init(4));
    deg2rad(x_init(5));
    deg2rad(x_final(4));
    deg2rad(x_final(5));
    deg2rad(gimbal_max);
    deg2rad(theta_max);
    deg2rad(gamma_gs);
}

void RocketLanding2D::systemFlowMap(
    const state_vector_ad_t &x,
    const input_vector_ad_t &u,
    state_vector_ad_t &f)
{
    typedef scalar_ad_t T;

    f(0) = x(2);
    f(1) = x(3);
    f(2) = 1. / T(m) * sin(x(4) + u(1)) * u(0);
    f(3) = 1. / T(m) * (cos(x(4) + u(1)) * u(0) - T(m) * T(g));
    f(4) = x(5);
    f(5) = 1. / T(I) * (-sin(u(1)) * u(0) * T(r_T));
}

void RocketLanding2D::initializeTrajectory(Eigen::MatrixXd &X,
                                           Eigen::MatrixXd &U)
{
    const size_t K = X.cols();

    X.setZero();
    U.setZero();
    for (size_t k = 0; k < K; k++)
    {
        const double alpha1 = double(K - k) / K;
        const double alpha2 = double(k) / K;

        X(0, k) = x_init(0) * alpha1 + x_final(0) * alpha2;
        X(1, k) = x_init(1) * alpha1 + x_final(1) * alpha2;
        X(2, k) = x_init(2) * alpha1 + x_final(2) * alpha2;
        X(3, k) = x_init(3) * alpha1 + x_final(3) * alpha2;
        X(4, k) = x_init(4) * alpha1 + x_final(4) * alpha2;
        X(5, k) = x_init(5) * alpha1 + x_final(5) * alpha2;

        U(0, k) = (T_max - T_min) / 2.;
        U(1, k) = 0;
    }
}

void RocketLanding2D::addApplicationConstraints(
    op::SecondOrderConeProgram &socp,
    Eigen::MatrixXd &X0,
    Eigen::MatrixXd &U0)
{
    const size_t K = X0.cols();

    auto var = [&](const string &name, const vector<size_t> &indices = {}) { return socp.get_variable(name, indices); };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };

    // initial state
    socp.addConstraint((-1.0) * var("X", {0, 0}) + param(x_init(0)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {1, 0}) + param(x_init(1)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {2, 0}) + param(x_init(2)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {3, 0}) + param(x_init(3)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {4, 0}) + param(x_init(4)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {5, 0}) + param(x_init(5)) == 0.0);
    // final state
    socp.addConstraint((-1.0) * var("X", {0, K - 1}) + param(x_final(0)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {1, K - 1}) + param(x_final(1)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {2, K - 1}) + param(x_final(2)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {3, K - 1}) + param(x_final(3)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {4, K - 1}) + param(x_final(4)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {5, K - 1}) + param(x_final(5)) == 0.0);

    for (size_t k = 0; k < K; ++k)
    {
        // glide slope
        socp.addConstraint(op::norm2({(1.0) * var("X", {0, k})}) <= (1.0 / tan(gamma_gs)) * var("X", {1, k}));

        // angle constraint
        socp.addConstraint((1.0) * var("X", {4, k}) + (theta_max) >= (0.0));
        socp.addConstraint((-1.0) * var("X", {4, k}) + (theta_max) >= (0.0));

        // throttle control constraints
        socp.addConstraint((1.0) * var("U", {0, k}) + (-T_min) >= (0.0));
        socp.addConstraint((-1.0) * var("U", {0, k}) + (T_max) >= (0.0));

        // gimbal control constraints
        socp.addConstraint((1.0) * var("U", {1, k}) + (gimbal_max) >= (0.0));
        socp.addConstraint((-1.0) * var("U", {1, k}) + (gimbal_max) >= (0.0));
    }
}

void RocketLanding2D::nondimensionalize()
{
    const double r_scale = x_init.segment(0, 2).norm();
    const double m_scale = m;

    r_T /= r_scale;
    g /= r_scale;
    I /= m_scale * r_scale * r_scale;
    m /= m_scale;
    T_min /= m_scale * r_scale;
    T_max /= m_scale * r_scale;

    x_init.segment(0, 4) /= r_scale;
    x_final.segment(0, 4) /= r_scale;
}

void RocketLanding2D::getStateWeightVector(state_vector_t &w)
{
    w.setOnes();
}

void RocketLanding2D::getInputWeightVector(input_vector_t &w)
{
    w.setOnes();
}

} // namespace rocket2d