#include "common.hpp"
#include "rocketLanding3d.hpp"

using std::string;
using std::vector;

namespace rocket3d
{

RocketLanding3D::RocketLanding3D()
{
    p.loadFromFile();
}

void RocketLanding3D::systemFlowMap(
    const state_vector_ad_t &x,
    const input_vector_ad_t &u,
    const param_vector_ad_t &p,
    state_vector_ad_t &f)
{
    typedef scalar_ad_t T;

    auto alpha_m_ = p(0);
    auto g_I_ = p.segment<3>(1);
    auto J_B_inv = p.segment<3>(4).asDiagonal().inverse();
    auto r_T_B_ = p.segment<3>(7);
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

void RocketLanding3D::getInitializedTrajectory(Eigen::MatrixXd &X,
                                               Eigen::MatrixXd &U)
{
    const size_t K = X.cols();

    for (size_t k = 0; k < K; k++)
    {
        const double alpha1 = double(K - k) / K;
        const double alpha2 = double(k) / K;

        // mass, position and linear velocity
        X(0, k) = alpha1 * p.x_init(0) + alpha2 * p.x_final(0);
        X.col(k).segment(1, 6) = alpha1 * p.x_init.segment(1, 6) + alpha2 * p.x_final.segment(1, 6);

        // do SLERP for quaternion
        Eigen::Quaterniond q0(p.x_init(7), p.x_init(8), p.x_init(9), p.x_init(10));
        Eigen::Quaterniond q1(p.x_final(7), p.x_final(8), p.x_final(9), p.x_final(10));
        Eigen::Quaterniond qs = q0.slerp(alpha2, q1);
        X.col(k).segment(7, 4) << qs.w(), qs.vec();

        // angular velocity
        X.col(k).segment(11, 3) = alpha1 * p.x_init.segment(11, 3) + alpha2 * p.x_final.segment(11, 3);

        // input
        U.setConstant((p.T_max - p.T_min) / 2.);
    }
}

void RocketLanding3D::addApplicationConstraints(
    op::SecondOrderConeProgram &socp,
    Eigen::MatrixXd &X0,
    Eigen::MatrixXd &U0)
{
    const size_t K = X0.cols();

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
    // socp.addConstraint((-1.0) * var("X", {0, K - 1}) + param(p.x_final(0)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {1, K - 1}) + param(p.x_final(1)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {2, K - 1}) + param(p.x_final(2)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {3, K - 1}) + param(p.x_final(3)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {4, K - 1}) + param(p.x_final(4)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {5, K - 1}) + param(p.x_final(5)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {6, K - 1}) + param(p.x_final(6)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {7, K - 1}) + param(p.x_final(7)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {8, K - 1}) + param(p.x_final(8)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {9, K - 1}) + param(p.x_final(9)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {10, K - 1}) + param(p.x_final(10)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {11, K - 1}) + param(p.x_final(11)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {12, K - 1}) + param(p.x_final(12)) == 0.0);
    // socp.addConstraint((-1.0) * var("X", {13, K - 1}) + param(p.x_final(13)) == 0.0);

    // Final Input
    socp.addConstraint((1.0) * var("U", {0, K - 1}) == (0.0));
    socp.addConstraint((1.0) * var("U", {1, K - 1}) == (0.0));

    // State Constraints:
    for (size_t k = 0; k < K; k++)
    {
        // Mass
        //     x(0) >= m_dry
        //     for all k
        socp.addConstraint((1.0) * var("X", {0, k}) + param_fn([this]() { return -p.x_final(0); }) >= (0.0));

        // Max Tilt Angle
        // norm2([x(8), x(9)]) <= sqrt((1 - cos_theta_max) / 2)
        socp.addConstraint(op::norm2({(1.0) * var("X", {8, k}),
                                      (1.0) * var("X", {9, k})}) <= param_fn([this]() { return sqrt((1.0 - cos(p.theta_max)) / 2.); }));

        // Glide Slope
        socp.addConstraint(
            op::norm2({(1.0) * var("X", {1, k}),
                       (1.0) * var("X", {2, k})}) <= param_fn([this]() { return tan(p.gamma_gs); }) * var("X", {3, k}));

        // Max Rotation Velocity
        socp.addConstraint(
            op::norm2({(1.0) * var("X", {11, k}),
                       (1.0) * var("X", {12, k}),
                       (1.0) * var("X", {13, k})}) <= param(p.w_B_max));
    }

    // Control Constraints
    for (size_t k = 0; k < K; k++)
    {
        // // Linearized Minimum Thrust
        op::AffineExpression lhs;
        for (size_t i = 0; i < INPUT_DIM_; i++)
        {
            lhs = lhs + param_fn([&U0, i, k]() { return (U0(i, k) / sqrt(U0(0, k) * U0(0, k) + U0(1, k) * U0(1, k) + U0(2, k) * U0(2, k))); }) * var("U", {i, k});
        }
        socp.addConstraint(lhs + param_fn([this]() { return -p.T_min; }) >= (0.0));

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

void RocketLanding3D::nondimensionalize()
{
    p.nondimensionalize();
}

void RocketLanding3D::redimensionalize()
{
    p.redimensionalize();
}

void RocketLanding3D::redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                Eigen::MatrixXd &U)
{
    p.redimensionalizeTrajectory(X, U);
}

void RocketLanding3D::getNewModelParameters(param_vector_t &param)
{
    param << p.alpha_m, p.g_I, p.J_B, p.r_T_B;
}

} // namespace rocket3d