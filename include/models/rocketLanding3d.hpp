#pragma once

#include <string>

#include <Eigen/Dense>

#include "systemModel.hpp"
#include "ecosWrapper.hpp"
#include "parameterServer.hpp"

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
        ParameterServer param(format("../include/models/config/{}.info", getModelName()));

        double I_sp;
        double m_init, m_dry;
        Eigen::Vector3d r_init, v_init, rpy_init, w_init;
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
        param.loadScalar("final_time_guess", final_time_guess);
        param.loadScalar("T_min", T_min);
        param.loadScalar("T_max", T_max);
        param.loadScalar("I_sp", I_sp);
        param.loadScalar("gimbal_max", gimbal_max);
        param.loadScalar("theta_max", theta_max);
        param.loadScalar("gamma_gs", gamma_gs);
        param.loadScalar("w_B_max", w_B_max);

        deg2rad(gimbal_max);
        deg2rad(theta_max);
        deg2rad(gamma_gs);
        deg2rad(w_init);
        deg2rad(w_B_max);
        deg2rad(rpy_init);

        alpha_m = 1. / (I_sp * fabs(g_I(2)));

        x_init << m_init, r_init, v_init, eulerToQuaternion(rpy_init), w_init;
        x_final << m_dry, r_final, v_final, 1., 0., 0., 0., 0, 0, 0;
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

        auto m = x(0);
        auto v_I = x.template segment<3>(4);
        auto q_B_I = x.template segment<4>(7);
        auto w_B = x.template segment<3>(11);

        f(0) = -alpha_m_ * u.norm();
        f.segment(1, 3) << v_I;
        f.segment(4, 3) << 1. / m * dirCosineMatrix<T>(q_B_I).transpose() * u + g_I_;
        f.segment(7, 4) << T(0.5) * omegaMatrix<T>(w_B) * q_B_I;
        f.segment(11, 3) << J_B_.inverse() * (r_T_B_.cross(u)) - w_B.cross(w_B);
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
        for (size_t i = 0; i < STATE_DIM_; i++)
        {
            socp.add_constraint((-1.0) * var("X", {i, 0}) + (x_init(i)) == 0.0);
        }

        // Final State
        // mass is free, quaternion
        // socp.add_constraint((-1.0) * var("X", {0, K - 1}) + (x_final(0)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {1, K - 1}) + (x_final(1)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {2, K - 1}) + (x_final(2)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {3, K - 1}) + (x_final(3)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {4, K - 1}) + (x_final(4)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {5, K - 1}) + (x_final(5)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {6, K - 1}) + (x_final(6)) == 0.0);
        // socp.add_constraint((-1.0) * var("X", {7, K - 1}) + (x_final(7)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {8, K - 1}) + (x_final(8)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {9, K - 1}) + (x_final(9)) == 0.0);
        // socp.add_constraint((-1.0) * var("X", {10, K - 1}) + (x_final(10)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {11, K - 1}) + (x_final(11)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {12, K - 1}) + (x_final(12)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {13, K - 1}) + (x_final(13)) == 0.0);

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

            // Maximum Thrust
            socp.add_constraint(
                op::norm2({(1.0) * var("U", {0, k}),
                           (1.0) * var("U", {1, k}),
                           (1.0) * var("U", {2, k})}) <= (T_max));

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

    state_vector_t getStateWeightVector() override
    {
        state_vector_t w;
        w.setOnes();
        return w;
    }

    input_vector_t getInputWeightVector() override
    {
        input_vector_t w;
        w.setOnes();
        w.head(3) *= T_max;
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

    state_vector_t x_init;
    state_vector_t x_final;

    double final_time_guess;

    double gimbal_max;
    double theta_max;
    double gamma_gs;
    double w_B_max;
};

} // namespace rocket3d