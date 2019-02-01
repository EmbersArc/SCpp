#pragma once

#include <string>

#include <Eigen/Dense>

#include "systemModel.hpp"
#include "ecosWrapper.hpp"

#include "rocketLanding2dDefinitions.hpp"

using std::string;

namespace rocket2d
{

class RocketLanding2D : public SystemModel<STATE_DIM_, INPUT_DIM_>
{
  public:
    RocketLanding2D();

    static string getModelName()
    {
        return "RocketLanding2D";
    }

    double getFinalTimeGuess() override
    {
        return final_time_guess;
    }

    void systemFlowMap(
        const state_vector_ad_t &x,
        const input_vector_ad_t &u,
        state_vector_ad_t &f) override;

    void initializeTrajectory(Eigen::MatrixXd &X,
                              Eigen::MatrixXd &U) override;

    void addApplicationConstraints(
        op::SecondOrderConeProgram &socp,
        Eigen::MatrixXd &X0,
        Eigen::MatrixXd &U0) override;

    void nondimensionalize() override;

    void redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                    Eigen::MatrixXd &U) override;

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

    state_vector_t x_init;
    state_vector_t x_final;
    double final_time_guess;
};

} // namespace rocket2d