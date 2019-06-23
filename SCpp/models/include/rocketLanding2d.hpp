#pragma once

#include <string>

#include "systemModel.hpp"
#include "ecosWrapper.hpp"

#include "rocketLanding2dDefinitions.hpp"

namespace rocket2d
{

/**
 * @brief A 2D rocket landing model.
 * 
 */
class RocketLanding2D : public SystemModel<STATE_DIM_, INPUT_DIM_, PARAM_DIM_>
{
  public:
    RocketLanding2D();

    static std::string getModelName()
    {
        return "RocketLanding2D";
    }

    void getTimeHorizon(double &T) override
    {
        T = final_time_guess;
    }

    void systemFlowMap(
        const state_vector_ad_t &x,
        const input_vector_ad_t &u,
        const param_vector_ad_t &par,
        state_vector_ad_t &f) override;

    void getInitializedTrajectory(Eigen::MatrixXd &X,
                                  Eigen::MatrixXd &U) override;

    void addApplicationConstraints(
        op::SecondOrderConeProgram &socp,
        Eigen::MatrixXd &X0,
        Eigen::MatrixXd &U0) override;

    void nondimensionalize() override;

    void redimensionalize() override;

    void nondimensionalizeTrajectory(Eigen::MatrixXd &X,
                                     Eigen::MatrixXd &U) override;

    void redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                    Eigen::MatrixXd &U) override;

    void getNewModelParameters(param_vector_t &param);

  private:
    double m;
    double g;
    double r_T;
    double I;
    double T_min;
    double T_max;
    double gimbal_max;
    double theta_max;
    double gamma_gs;

    double m_scale, r_scale;

    state_vector_t x_init;
    state_vector_t x_final;
    double final_time_guess;
};

} // namespace rocket2d