#pragma once

#include <string>

#include <Eigen/Dense>

#include "systemModel.hpp"
#include "ecosWrapper.hpp"
#include "parameters.hpp"

#include "rocketLanding3dDefinitions.hpp"

using std::string;

namespace rocket3d
{

class RocketLanding3D : public SystemModel<RocketLanding3D, STATE_DIM_, INPUT_DIM_>
{
  public:
    typedef SystemModel<RocketLanding3D, STATE_DIM_, INPUT_DIM_> BASE;

    RocketLanding3D()
    {
        string configFilePath = format("../include/models/config/{}.info", getModelName());
        loadMatrix(configFilePath, "g_I", g_I);
        loadMatrix(configFilePath, "J_B", J_B);
        loadMatrix(configFilePath, "r_T_B", r_T_B);

        loadMatrix(configFilePath, "x_init", x_init);
        loadMatrix(configFilePath, "x_final", x_final);

        loadScalar(configFilePath, "constrain_initial_orientation", constrain_initial_orientation);

        if (!constrain_initial_orientation)
        {
            x_init.segment(8, 3) << 1., 0., 0., 0.;
        }

        loadScalar(configFilePath, "T_min", T_min);
        loadScalar(configFilePath, "T_max", T_max);
        loadScalar(configFilePath, "alpha_m", alpha_m);

        double delta_max;
        double theta_max;
        double gamma_gs;

        loadScalar(configFilePath, "delta_max", delta_max);
        loadScalar(configFilePath, "theta_max", theta_max);
        loadScalar(configFilePath, "gamma_gs", gamma_gs);
        loadScalar(configFilePath, "w_B_max", w_B_max);

        deg2rad(delta_max);
        deg2rad(theta_max);
        deg2rad(gamma_gs);
        deg2rad(w_B_max);

        tan_delta_max = tan(delta_max);
        cos_theta_max = cos(theta_max);
        tan_gamma_gs = tan(gamma_gs);
    }

    static string getModelName()
    {
        return "RocketLanding3D";
    }

    static double getFinalTimeGuess()
    {
        return 8.;
    }

    template <typename T>
    void systemFlowMap(
        const Eigen::Matrix<T, BASE::state_dim_, 1> &x,
        const Eigen::Matrix<T, BASE::input_dim_, 1> &u,
        Eigen::Matrix<T, BASE::state_dim_, 1> &f)
    {
        auto alpha_m_ = T(alpha_m);
        auto g_I_ = g_I.cast<T>();
        auto J_B_ = J_B.cast<T>().asDiagonal();
        auto r_T_B_ = r_T_B.cast<T>();

        f(0) = -alpha_m_ * u.norm();
        f.segment(1, 3) << x.segment(4, 3);
        f.segment(4, 3) << 1. / x(0) * dirCosineMatrix<T>(x.segment(7, 4)).transpose() * u + g_I_;
        f.segment(7, 4) << T(0.5) * omegaMatrix<T>(x.segment(11, 3)) * x.segment(7, 4);
        f.segment(11, 3) << J_B_.inverse() * (skew<T>(r_T_B_) * u - skew<T>(x.segment(11, 3)) * J_B_ * x.segment(11, 3));
    }

    void initializeTrajectory(Eigen::MatrixXd &X,
                              Eigen::MatrixXd &U) override
    {
        //    Nondimensionalize();

        const size_t K = X.cols();

        for (size_t k = 0; k < K; k++)
        {
            const double alpha1 = double(K - k) / K;
            const double alpha2 = double(k) / K;

            // mass, position and linear velocity
            X(0, k) = alpha1 * x_init(0) + alpha2 * x_final(0);
            X.col(k).segment(1, 6) = alpha1 * x_init.segment(1, 6) + alpha2 * x_final.segment(1, 6);

            // do SLERP for quaternion
            Eigen::Quaterniond q0, q1;
            q0.w() = x_init(7);
            q0.vec() = x_init.segment(8, 3);
            q1.w() = x_final(0);
            q1.vec() << x_final.segment(8, 3);
            Eigen::Quaterniond slerpQuaternion = q0.slerp(alpha2, q1);
            X.col(k).segment(7, 4) << slerpQuaternion.w(), slerpQuaternion.vec();

            // angular velocity
            X.col(k).segment(11, 3) = alpha1 * x_init.segment(11, 3) + alpha2 * x_final.segment(11, 3);

            // input
            U.col(k) = (alpha1 * x_init(0) + alpha2 * x_final(0)) * -g_I;
        }
    }

    void addApplicationConstraints(
        optimization_problem::SecondOrderConeProgram &socp,
        Eigen::MatrixXd &X0,
        Eigen::MatrixXd &U0) override
    {
        const size_t K = X0.cols();

        auto var = [&](const string &name, const vector<size_t> &indices) { return socp.get_variable(name, indices); };
        //    auto param = [](double &param_value){ return optimization_problem::Parameter(&param_value); };
        //    auto param_fn = [](std::function<double()> callback){ return optimization_problem::Parameter(callback); };

        // Initial state
        socp.add_constraint((-1.0) * var("X", {0, 0}) + (x_init(0)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {1, 0}) + (x_init(1)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {2, 0}) + (x_init(2)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {3, 0}) + (x_init(3)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {4, 0}) + (x_init(4)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {5, 0}) + (x_init(5)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {6, 0}) + (x_init(6)) == 0.0);
        if (constrain_initial_orientation)
        {
            socp.add_constraint((-1.0) * var("X", {7, 0}) + (x_init(7)) == 0.0);
            socp.add_constraint((-1.0) * var("X", {8, 0}) + (x_init(8)) == 0.0);
            socp.add_constraint((-1.0) * var("X", {9, 0}) + (x_init(9)) == 0.0);
            socp.add_constraint((-1.0) * var("X", {10, 0}) + (x_init(10)) == 0.0);
        }
        socp.add_constraint((-1.0) * var("X", {11, 0}) + (x_init(11)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {12, 0}) + (x_init(12)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {13, 0}) + (x_init(13)) == 0.0);

        // Final State (mass is free)
        for (size_t i = 1; i < STATE_DIM_; i++)
        {
            socp.add_constraint((-1.0) * var("X", {i, K - 1}) + (x_final(i)) == 0.0);
        }
        socp.add_constraint((1.0) * var("U", {0, K - 1}) == (0.0));
        socp.add_constraint((1.0) * var("U", {1, K - 1}) == (0.0));

        // State Constraints:
        for (size_t k = 0; k < K; k++)
        {
            // Mass
            //     x(0) >= m_dry
            //     for all k
            socp.add_constraint((1.0) * var("X", {0, k}) + (-x_final(0)) >= (0.0));

            // Max Tilt Angle
            //
            // norm2([x(8), x(9)]) <= c
            // with c := sqrt((1 - cos_theta_max) / 2)
            const double c = sqrt((1.0 - cos_theta_max) / 2.);
            socp.add_constraint(optimization_problem::norm2({(1.0) * var("X", {8, k}),
                                                             (1.0) * var("X", {9, k})}) <= (c));

            // Glide Slope
            socp.add_constraint(
                optimization_problem::norm2({(1.0) * var("X", {1, k}),
                                             (1.0) * var("X", {2, k})}) <= (1.0 / tan_gamma_gs) * var("X", {3, k}));

            // Max Rotation Velocity
            socp.add_constraint(
                optimization_problem::norm2({(1.0) * var("X", {11, k}),
                                             (1.0) * var("X", {12, k}),
                                             (1.0) * var("X", {13, k})}) <= (w_B_max));
        }

        // Control Constraints
        for (size_t k = 0; k < K; k++)
        {
            // Linearized Minimum Thrust
            // optimization_problem::AffineExpression lhs;
            // for (size_t i = 0; i < n_inputs; i++)
            // {
            //     lhs = lhs + param_fn([&U0, i, k]() { return (U0(i, k) / sqrt(U0(0, k) * U0(0, k) + U0(1, k) * U0(1, k) + U0(2, k) * U0(2, k))); }) * var("U", {i, k});
            // }
            // socp.add_constraint(lhs + (-T_min) >= (0.0));

            // Simplified Minimum Thrust
            socp.add_constraint((1.0) * var("U", {2, k}) + (-T_min) >= (0.0));

            // Maximum Thrust
            socp.add_constraint(
                optimization_problem::norm2({(1.0) * var("U", {0, k}),
                                             (1.0) * var("U", {1, k}),
                                             (1.0) * var("U", {2, k})}) <= (T_max));

            // Maximum Gimbal Angle
            socp.add_constraint(
                optimization_problem::norm2({(1.0) * var("U", {0, k}),
                                             (1.0) * var("U", {1, k})}) <= (tan_delta_max)*var("U", {2, k}));
        }
    }

    // void Nondimensionalize()
    // {
    //     double r_scale = r_I_init.norm();
    //     double m_scale = m_wet;

    //     alpha_m *= r_scale;
    //     r_T_B /= r_scale;
    //     g_I /= r_scale;
    //     J_B /= (m_scale * r_scale * r_scale);

    //     m_wet /= m_scale;
    //     r_I_init /= r_scale;
    //     v_I_init /= r_scale;

    //     m_dry /= m_scale;
    //     r_I_final /= r_scale;
    //     v_I_final /= r_scale;

    //     T_max /= m_scale * r_scale;
    //     T_min /= m_scale * r_scale;
    // }

  private:
    Eigen::Vector3d g_I;
    Eigen::Vector3d J_B;
    Eigen::Vector3d r_T_B;
    double alpha_m;
    double T_min;
    double T_max;

    state_vector_t x_init;
    state_vector_t x_final;

    bool constrain_initial_orientation;

    double tan_delta_max;
    double cos_theta_max;
    double tan_gamma_gs;
    double w_B_max;
};

} // namespace rocket3d