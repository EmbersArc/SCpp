#pragma once

#include <string>

#include <Eigen/Dense>

#include "systemModel.hpp"
#include "ecosWrapper.hpp"

#include "rocketLanding3dDefinitions.hpp"

using std::string;

namespace rocket3d
{

class RocketLanding3D : public SystemModel<STATE_DIM_, INPUT_DIM_>
{
  public:
    typedef SystemModel<STATE_DIM_, INPUT_DIM_> BASE;

    RocketLanding3D();

    static string getModelName()
    {
        return "RocketLanding3D";
    }

    double getFinalTimeGuess() const
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

    void nondimensionalize();

    void redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                    Eigen::MatrixXd &U);

    void getStateWeightVector(state_vector_t &w) override;
    void getInputWeightVector(input_vector_t &w) override;

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