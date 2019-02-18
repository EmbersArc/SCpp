#pragma once

#include <string>
#include <random>

#include "systemModel.hpp"
#include "ecosWrapper.hpp"
#include "parameterServer.hpp"

#include "rocketLanding3dDefinitions.hpp"

namespace rocket3d
{

/**
 * @brief A 3D rocket landing model.
 * 
 */
class RocketLanding3D : public SystemModel<STATE_DIM_, INPUT_DIM_, PARAM_DIM_>
{
  public:
    RocketLanding3D();

    static std::string getModelName()
    {
        return "RocketLanding3D";
    }

    void getFinalTimeGuess(double &sigma) override
    {
        sigma = p.final_time_guess;
    }

    void systemFlowMap(
        const state_vector_ad_t &x,
        const input_vector_ad_t &u,
        const param_vector_ad_t &p,
        state_vector_ad_t &f) override;

    void getInitializedTrajectory(Eigen::MatrixXd &X,
                                  Eigen::MatrixXd &U) override;

    void addApplicationConstraints(
        op::SecondOrderConeProgram &socp,
        Eigen::MatrixXd &X0,
        Eigen::MatrixXd &U0) override;

    void nondimensionalize() override;

    void redimensionalize() override;

    void redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                    Eigen::MatrixXd &U) override;

    void getNewModelParameters(param_vector_t &param);

    /**
     * @brief Varies the initial state randomly.
     * 
     */
    void randomizeInitialState();

    struct RocketLanding3dParameters
    {
        Eigen::Vector3d g_I;
        Eigen::Vector3d J_B;
        Eigen::Vector3d r_T_B;
        double alpha_m;
        double T_min;
        double T_max;

        double gimbal_max;
        double theta_max;
        double gamma_gs;
        double w_B_max;

        Eigen::Vector3d rpy_init;

        state_vector_t x_init;
        state_vector_t x_final;
        double final_time_guess;

        double m_scale, r_scale;

        void randomizeInitialState()
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
            x_init.segment(7, 4) << eulerToQuaternion(euler);
        }

        void loadFromFile()
        {
            ParameterServer param(fmt::format("../models/config/{}.info", getModelName()));

            bool randomInitialState;
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
            param.loadScalar("final_time_guess", final_time_guess);
            param.loadScalar("T_min", T_min);
            param.loadScalar("T_max", T_max);
            param.loadScalar("I_sp", I_sp);
            param.loadScalar("gimbal_max", gimbal_max);
            param.loadScalar("theta_max", theta_max);
            param.loadScalar("gamma_gs", gamma_gs);
            param.loadScalar("w_B_max", w_B_max);
            param.loadScalar("randomInitialState", randomInitialState);

            deg2rad(gimbal_max);
            deg2rad(theta_max);
            deg2rad(gamma_gs);
            deg2rad(w_init);
            deg2rad(w_B_max);
            deg2rad(rpy_init);

            alpha_m = 1. / (I_sp * fabs(g_I(2)));

            x_init << m_init, r_init, v_init, eulerToQuaternion(rpy_init), w_init;
            if (randomInitialState)
            {
                randomizeInitialState();
            }
            x_final << m_dry, r_final, v_final, 1., 0., 0., 0., 0, 0, 0;
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

        void redimensionalize()
        {
            alpha_m /= r_scale;
            r_T_B *= r_scale;
            g_I *= r_scale;
            J_B *= (m_scale * r_scale * r_scale);

            x_init(0) *= m_scale;
            x_init.segment(1, 3) *= r_scale;
            x_init.segment(4, 3) *= r_scale;

            x_final(0) *= m_scale;
            x_final.segment(1, 3) *= r_scale;
            x_final.segment(4, 3) *= r_scale;

            T_min *= m_scale * r_scale;
            T_max *= m_scale * r_scale;
        }

        void redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                        Eigen::MatrixXd &U)
        {
            const size_t K = X.cols();

            X.row(0) *= m_scale;
            X.block(1, 0, 6, K) *= r_scale;

            U *= m_scale * r_scale;
        }
    } p;

  private:
};

} // namespace rocket3d