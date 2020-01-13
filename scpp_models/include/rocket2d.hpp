#pragma once

#include <string>
#include <random>

#include "systemModel.hpp"
#include "socpSolver.hpp"
#include "parameterServer.hpp"

#include "rocket2dDefinitions.hpp"

namespace scpp::models
{

/**
 * @brief A 2D rocket model.
 * 
 */
class Rocket2d : public SystemModel<Rocket2d, STATE_DIM_, INPUT_DIM_, PARAM_DIM_>
{
public:
    Rocket2d();

    inline static const std::string modelName = "Rocket2D";

    void systemFlowMap(
        const state_vector_ad_t &x,
        const input_vector_ad_t &u,
        const param_vector_ad_t &par,
        state_vector_ad_t &f) override;

    void getOperatingPoint(state_vector_t &x, input_vector_t &u) override;

    void addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                   state_vector_v_t &X0,
                                   input_vector_v_t &U0) override;

    void nondimensionalize() override;

    void redimensionalize() override;

    void nondimensionalizeTrajectory(trajectory_data_t &td) override;

    void redimensionalizeTrajectory(trajectory_data_t &td) override;

    void getNewModelParameters(param_vector_t &param) override;

    void getInitializedTrajectory(trajectory_data_t &td) override;

    void loadParameters();

    struct Parameters
    {
        double m;
        double J_B;
        Eigen::Vector2d g_I;
        Eigen::Vector2d r_T_B;
        double T_min;
        double T_max;

        double gamma_gs;
        double gimbal_max;
        double theta_max;
        double v_I_max;
        double w_B_max;

        double eta_init;
        double eta_final;

        state_vector_t x_init;
        state_vector_t x_final;
        double final_time;

        bool constrain_initial_final;
        bool add_slack_variables;

        void loadFromFile(const std::string &path);

        void nondimensionalize();

        void redimensionalize();
    } p;
};

} // namespace scpp::models