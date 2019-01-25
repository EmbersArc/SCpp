#pragma once

#include <string>

#include <Eigen/Dense>

#include "systemModel.hpp"
#include "ecosWrapper.hpp"
#include "parameters.hpp"

#include "rocketLanding2dDefinitions.hpp"

using std::string;

namespace rocket2d
{

class RocketLanding2D : public SystemModel<RocketLanding2D, STATE_DIM_, INPUT_DIM_>
{
  public:
    typedef SystemModel<RocketLanding2D, STATE_DIM_, INPUT_DIM_> BASE;

    RocketLanding2D()
    {
        string configFilePath = format("../include/models/config/{}.info", getModelName());

        loadScalar(configFilePath, "m", m);
        loadScalar(configFilePath, "g", g);
        loadScalar(configFilePath, "r_T", r_T);
        loadScalar(configFilePath, "I", I);
        loadScalar(configFilePath, "T_min", T_min);
        loadScalar(configFilePath, "T_max", T_max);
        loadScalar(configFilePath, "gimbal_max", gimbal_max);

        loadMatrix(configFilePath, "x_init", x_init);
        loadMatrix(configFilePath, "x_final", x_final);

        nondimensionalize();
    }

    static string getModelName()
    {
        return "RocketLanding2D";
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
        f(0) = x(2);
        f(1) = x(3);
        f(2) = 1. / T(m) * sin(x(4) + u(1)) * u(0);
        f(3) = 1. / T(m) * (cos(x(4) + u(1)) * u(0) - T(m) * T(g));
        f(4) = x(5);
        f(5) = 1. / I * (-sin(u(1)) * u(0) * T(r_T));
    }

    void initializeTrajectory(Eigen::MatrixXd &X,
                              Eigen::MatrixXd &U) override
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

    void addApplicationConstraints(
        op::SecondOrderConeProgram &socp,
        Eigen::MatrixXd &X0,
        Eigen::MatrixXd &U0) override
    {
        const size_t K = X0.cols();

        auto var = [&](const string &name, const vector<size_t> &indices) { return socp.get_variable(name, indices); };
        auto param = [](double &param_value) { return op::Parameter(&param_value); };

        // initial state
        socp.add_constraint((-1.0) * var("X", {0, 0}) + param(x_init(0)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {1, 0}) + param(x_init(1)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {2, 0}) + param(x_init(2)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {3, 0}) + param(x_init(3)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {4, 0}) + param(x_init(4)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {5, 0}) + param(x_init(5)) == 0.0);
        // final state
        socp.add_constraint((-1.0) * var("X", {0, K - 1}) + param(x_final(0)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {1, K - 1}) + param(x_final(1)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {2, K - 1}) + param(x_final(2)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {3, K - 1}) + param(x_final(3)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {4, K - 1}) + param(x_final(4)) == 0.0);
        socp.add_constraint((-1.0) * var("X", {5, K - 1}) + param(x_final(5)) == 0.0);

        // glide slope cone
        // TODO

        for (size_t k = 0; k < K; ++k)
        {
            // throttle control constraints
            socp.add_constraint((1.0) * var("U", {0, k}) + (-T_min) >= (0.0));
            socp.add_constraint((-1.0) * var("U", {0, k}) + (T_max) >= (0.0));

            // gimbal control constraints
            socp.add_constraint((1.0) * var("U", {1, k}) + (gimbal_max) >= (0.0));
            socp.add_constraint((-1.0) * var("U", {1, k}) + (gimbal_max) >= (0.0));
        }
    }

    void nondimensionalize()
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
        return w;
    }

  private:
    double m;
    double g;
    double r_T;
    double I;
    double T_min;
    double T_max;
    double gimbal_max;

    state_vector_t x_init;
    state_vector_t x_final;
};

} // namespace rocket2d