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
        loadScalar(configFilePath, "final_time_guess", final_time_guess);

        if (!constrain_initial_orientation)
        {
            x_init.segment(7, 4) << 1., 0., 0., 0.;
        }

        loadScalar(configFilePath, "T_min", T_min);
        loadScalar(configFilePath, "T_max", T_max);
        loadScalar(configFilePath, "t_max", t_max);

        double I_sp;
        loadScalar(configFilePath, "I_sp", I_sp);
        alpha_m = 1. / (I_sp * fabs(g_I(2)));

        loadScalar(configFilePath, "gimbal_max", gimbal_max);
        loadScalar(configFilePath, "theta_max", theta_max);
        loadScalar(configFilePath, "gamma_gs", gamma_gs);
        loadScalar(configFilePath, "w_B_max", w_B_max);

        deg2rad(gimbal_max);
        deg2rad(theta_max);
        deg2rad(gamma_gs);
        deg2rad(w_B_max);
    }

    static string getModelName()
    {
        return "RocketLanding3D";
    }

    double getFinalTimeGuess() const
    {
        return final_time_guess;
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

        auto T_B = u.template head<3>();
        // Eigen::Matrix<T, 3, 1> T_r;
        // T_r << T(0.), T(0.), u(3);

        auto m = x(0);
        auto v_I = x.template segment<3>(4);
        auto q_B_I = x.template segment<4>(7);
        auto w_B = x.template segment<3>(11);

        f(0) = -alpha_m_ * T_B.norm();
        f.segment(1, 3) << v_I;
        f.segment(4, 3) << 1. / m * dirCosineMatrix<T>(q_B_I).transpose() * T_B + g_I_;
        f.segment(7, 4) << T(0.5) * omegaMatrix<T>(w_B) * q_B_I;
        f.segment(11, 3) << J_B_.inverse() * (r_T_B_.cross(T_B)) - w_B.cross(w_B);
    }

    void initializeTrajectory(Eigen::MatrixXd &X,
                              Eigen::MatrixXd &U) override
    {
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
            U.col(k).head(3) = (alpha1 * x_init(0) + alpha2 * x_final(0)) * -g_I;
            // U.col(k)(3) = 0.;
        }
    }

    void addApplicationConstraints(
        op::SecondOrderConeProgram &socp,
        Eigen::MatrixXd &X0,
        Eigen::MatrixXd &U0) override
    {
        const size_t K = X0.cols();

        auto var = [&](const string &name, const vector<size_t> &indices = {}) { return socp.get_variable(name, indices); };
        // auto param = [](double &param_value){ return op::Parameter(&param_value); };
        auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

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
            // norm2([x(8), x(9)]) <= sqrt((1 - cos_theta_max) / 2)
            socp.add_constraint(op::norm2({(1.0) * var("X", {8, k}),
                                           (1.0) * var("X", {9, k})}) <= sqrt((1.0 - cos(theta_max)) / 2.));

            // Glide Slope
            socp.add_constraint(
                op::norm2({(1.0) * var("X", {1, k}),
                           (1.0) * var("X", {2, k})}) <= (1.0 / tan(gamma_gs)) * var("X", {3, k}));

            // Max Rotation Velocity
            socp.add_constraint(
                op::norm2({(1.0) * var("X", {11, k}),
                           (1.0) * var("X", {12, k}),
                           (1.0) * var("X", {13, k})}) <= (w_B_max));
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
            socp.add_constraint(lhs + (-T_min) >= (0.0));

            // // Simplified Minimum Thrust
            // socp.add_constraint((1.0) * var("U", {2, k}) + (-T_min) >= (0.0));

            // Maximum Thrust
            socp.add_constraint(
                op::norm2({(1.0) * var("U", {0, k}),
                           (1.0) * var("U", {1, k}),
                           (1.0) * var("U", {2, k})}) <= (T_max));

            // Maximum Roll Torque
            // socp.add_constraint((-1.0) * var("U", {3, k}) + (t_max) >= (0.0));

            // Maximum Gimbal Angle
            socp.add_constraint(
                op::norm2({(1.0) * var("U", {0, k}),
                           (1.0) * var("U", {1, k})}) <= tan(gimbal_max) * var("U", {2, k}));
        }
    }

    void nondimensionalize()
    {
        m_scale = x_init(0);
        r_scale = x_init.segment(1, 3).norm();

        alpha_m *= r_scale;
        r_T_B /= r_scale;
        g_I /= r_scale;
        J_B /= (m_scale * r_scale * r_scale);

        x_init(0) /= m_scale;
        x_init.segment(1, 3) /= r_scale;
        x_init.segment(4, 3) /= r_scale;

        x_final(0) /= m_scale;
        x_final.segment(1, 3) /= r_scale;
        x_final.segment(4, 3) /= r_scale;

        T_min /= m_scale * r_scale;
        T_max /= m_scale * r_scale;

        t_max /= m_scale * r_scale * r_scale;
    }

    void redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                    Eigen::MatrixXd &U)
    {
        const size_t K = X.cols();

        X.row(0) *= m_scale;
        X.block(1, 0, 6, K) *= r_scale;

        U.topRows(3) *= m_scale * r_scale;
        U.bottomRows(1) *= m_scale * r_scale * r_scale;
    }

    state_vector_t getStateWeightVector()
    {
        state_vector_t w;
        w.setOnes();
        return w;
    }

    input_vector_t getInputWeightVector()
    {
        input_vector_t w;
        w.setOnes();
        w.head(3) *= T_max;
        w(3) *= t_max;
        return w;
    }

  private:
    double m_scale = 1.;
    double r_scale = 1.;

    Eigen::Vector3d g_I;
    Eigen::Vector3d J_B;
    Eigen::Vector3d r_T_B;
    double alpha_m;
    double T_min;
    double T_max;
    double t_max;

    state_vector_t x_init;
    state_vector_t x_final;

    bool constrain_initial_orientation;
    double final_time_guess;

    double gimbal_max;
    double theta_max;
    double gamma_gs;
    double w_B_max;
};

} // namespace rocket3d