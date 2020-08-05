#pragma once

#include <string>
#include <random>

#include "systemModel.hpp"
#include "parameterServer.hpp"
#include "common.hpp"

#include "rocketQuatDefinitions.hpp"

namespace scpp::models
{

/**
 * @brief A 3D rocket landing model.
 * 
 */
class RocketQuat : public SystemModel<RocketQuat, STATE_DIM, INPUT_DIM, PARAM_DIM>
{
public:
    RocketQuat() = default;

    inline static const std::string modelName = "RocketQuat";

    void systemFlowMap(
        const state_vector_ad_t &x,
        const input_vector_ad_t &u,
        const param_vector_ad_t &par,
        state_vector_ad_t &f) override;

    void getInitializedTrajectory(trajectory_data_t &td) override;

    void addApplicationConstraints(cvx::OptimizationProblem &socp,
                                   state_vector_v_t &X0,
                                   input_vector_v_t &U0) override;

    void nondimensionalize() override;

    void redimensionalize() override;

    void nondimensionalizeTrajectory(trajectory_data_t &td) override;

    void redimensionalizeTrajectory(trajectory_data_t &td) override;

    void getNewModelParameters(param_vector_t &param) override;

    void loadParameters();

    struct Parameters
    {
        bool exact_minimum_thrust;
        bool enable_roll_control;

        Eigen::Vector3d g_I;
        Eigen::Vector3d J_B;
        Eigen::Vector3d r_T_B;
        double alpha_m;
        double T_min;
        double T_max;
        double t_max;

        double gimbal_max;
        double theta_max;
        double gamma_gs;
        double w_B_max;

        state_vector_t x_init;
        state_vector_t x_final;
        double final_time;

        double m_scale, r_scale;

        void randomizeInitialState();

        void loadFromFile(const std::string &path);

        void nondimensionalize();

        void redimensionalize();

        void nondimensionalizeTrajectory(state_vector_v_t &X, input_vector_v_t &U) const;

        void redimensionalizeTrajectory(state_vector_v_t &X, input_vector_v_t &U) const;
    } p;

    struct DynamicParameters
    {
        double tilt_const;
        double gs_const;
        double gimbal_const;
        Eigen::MatrixXd thrust_const;
        input_vector_v_t *U0_ptr;
    } p_dyn;

private:
    void updateProblemParameters();
};

} // namespace scpp::models