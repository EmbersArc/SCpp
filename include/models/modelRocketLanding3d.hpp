#pragma once

#include <Eigen/Dense>

#include "modelRocketLanding3dDefinitions.hpp"

#include "systemModel.hpp"
#include "ecosWrapper.hpp"
#include "constants.hpp"

namespace rocket3d
{

template <typename T>
Eigen::Matrix<T, 3, 3> skew(const Eigen::Matrix<T, 3, 1> &v)
{
    Eigen::Matrix<T, 3, 3> skewMatrix;
    skewMatrix << T(0), -v(2), v(1),
        v(2), T(0), -v(0),
        -v(1), v(0), T(0);
    return skewMatrix;
}

template <typename T>
Eigen::Matrix<T, 3, 3> dirCosineMatrix(const Eigen::Matrix<T, 4, 1> &q)
{
    Eigen::Matrix<T, 3, 3> dirCosineMatrix;
    dirCosineMatrix << 1 - 2 * (q(2) * q(2) + q(3) * q(3)), 2 * (q(1) * q(2) + q(0) * q(3)), 2 * (q(1) * q(3) - q(0) * q(2)),
        2 * (q(1) * q(2) - q(0) * q(3)), 1 - 2 * (q(1) * q(1) + q(3) * q(3)), 2 * (q(2) * q(3) + q(0) * q(1)),
        2 * (q(1) * q(3) + q(0) * q(2)), 2 * (q(2) * q(3) - q(0) * q(1)), 1 - 2 * (q(1) * q(1) + q(2) * q(2));

    return dirCosineMatrix;
}

template <typename T>
Eigen::Matrix<T, 4, 4> omegaMatrix(const Eigen::Matrix<T, 3, 1> &w)
{
    Eigen::Matrix<T, 4, 4> omegaMatrix;
    omegaMatrix << T(0), -w(0), -w(1), -w(2),
        w(0), T(0), w(2), -w(1),
        w(1), -w(2), T(0), w(0),
        w(2), w(1), -w(0), T(0);

    return omegaMatrix;
}

class ModelRocketLanding3D : public SystemModel<ModelRocketLanding3D, rocket3d::STATE_DIM_, rocket3d::INPUT_DIM_, K>
{
  public:
    typedef SystemModel<ModelRocketLanding3D, rocket3d::STATE_DIM_, rocket3d::INPUT_DIM_, K> BASE;

    ModelRocketLanding3D(){};

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

    void initializeTrajectory(state_trajectory_matrix_t &X,
                              input_trajectory_matrix_t &U) override
    {
        //    Nondimensionalize();

        if (!constrain_initial_orientation)
        {
            q_B_I_init = {1., 0., 0., 0.};
        }

        state_vector_t x_init;
        x_init << m_wet, r_I_init, v_I_init, q_B_I_init, w_B_init;
        state_vector_t x_final;
        x_final << m_dry, r_I_final, v_I_final, q_B_I_final, w_B_final;

        for (int k = 0; k < K; k++)
        {
            const double alpha1 = double(K - k) / K;
            const double alpha2 = double(k) / K;

            // mass, position and linear velocity
            X(0, k) = alpha1 * x_init(0) + alpha2 * x_final(0);
            X.col(k).segment(1, 6) = alpha1 * x_init.segment(1, 6) + alpha2 * x_final.segment(1, 6);

            // do SLERP for quaternion
            Eigen::Quaterniond q0, q1;
            q0.w() = q_B_I_init(0);
            q0.vec() = q_B_I_init.tail<3>();
            q1.w() = q_B_I_final(0);
            q1.vec() << q_B_I_final.tail<3>();
            Eigen::Quaterniond slerpQuaternion = q0.slerp(alpha2, q1);
            X.col(k).segment(7, 4) << slerpQuaternion.w(), slerpQuaternion.vec();

            // angular velocity
            X.col(k).segment(11, 3) = alpha1 * x_init.segment(11, 3) + alpha2 * x_final.segment(11, 3);

            // input
            U.col(k) = (alpha1 * m_wet + alpha2 * m_dry) * -g_I;
        }
    }

    void addApplicationConstraints(
        optimization_problem::SecondOrderConeProgram &socp,
        state_trajectory_matrix_t &X0,
        input_trajectory_matrix_t &U0) override
    {
        auto var = [&](const string &name, const vector<size_t> &indices) { return socp.get_variable(name, indices); };
        //    auto param = [](double &param_value){ return optimization_problem::Parameter(&param_value); };
        //    auto param_fn = [](std::function<double()> callback){ return optimization_problem::Parameter(callback); };

        state_vector_t x_init;
        x_init << m_wet, r_I_init, v_I_init, q_B_I_init, w_B_init;
        state_vector_t x_final;
        x_final << m_dry, r_I_final, v_I_final, q_B_I_final, w_B_final;

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
        for (size_t i = 1; i < rocket3d::STATE_DIM_; i++)
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
            socp.add_constraint((1.0) * var("X", {0, k}) + (-m_dry) >= (0.0));

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

  private:
    string modelName_;

    Eigen::Vector3d g_I = {0, 0, -1};
    Eigen::Vector3d J_B = {0.168 * 2e-2, 0.168 * 1, 0.168 * 1};
    Eigen::Vector3d r_T_B = {0, 0, -0.25};
    double I_sp = 30.;
    double alpha_m = 1. / (I_sp * abs(g_I(2)));
    double T_min = 1.5;
    double T_max = 6.5;

    //initial state
    double m_wet = 2.;
    Eigen::Vector3d r_I_init = {0, 4, 4};
    Eigen::Vector3d v_I_init = {1, 0, 0};
    Eigen::Vector4d q_B_I_init = {1., 0., 0., 0.};
    bool constrain_initial_orientation = false;
    Eigen::Vector3d w_B_init = {0., 0., 0.};

    //final state
    double m_dry = 1.;
    Eigen::Vector3d r_I_final = {0., 0., 0.};
    Eigen::Vector3d v_I_final = {0, 0., -1e-1};
    Eigen::Vector4d q_B_I_final = {1., 0., 0., 0.};
    Eigen::Vector3d w_B_final = {0., 0., 0.};

    const double tan_delta_max = tan(20. / 180. * M_PI);
    const double cos_theta_max = cos(90. / 180. * M_PI);
    const double tan_gamma_gs = tan(75. / 180. * M_PI);
    const double w_B_max = 30. / 180. * M_PI;
};

} // namespace rocket3d