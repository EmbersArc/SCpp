#pragma once

#include <string>
#include <random>

#include "systemModel.hpp"
#include "ecosWrapper.hpp"
#include "parameterServer.hpp"

#include "rocketHoverDefinitions.hpp"

namespace rocketHover
{

/**
 * @brief A 3D rocket model.
 * 
 */
class RocketHover : public SystemModel<RocketHover, STATE_DIM_, INPUT_DIM_, PARAM_DIM_>
{
public:
    RocketHover();

    inline static const std::string modelName = "RocketHover";

    void systemFlowMap(
        const state_vector_ad_t &x,
        const input_vector_ad_t &u,
        const param_vector_ad_t &par,
        state_vector_ad_t &f) override;

    void getOperatingPoint(state_vector_t &x, input_vector_t &u) override;

    void addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                   Eigen::MatrixXd &X0,
                                   Eigen::MatrixXd &U0) override;

    void nondimensionalize() override;

    void redimensionalize() override;

    void getInitializedTrajectory(Eigen::MatrixXd &X,
                                  Eigen::MatrixXd &U) override;

    void loadParameters();

    struct Parameters
    {
        Eigen::Vector3d g_I;
        Eigen::Vector3d J_B;
        Eigen::Vector3d r_T_B;
        double m;
        double T_min;
        double T_max;

        double gimbal_max;
        double theta_max;
        double v_I_max;
        double w_B_max;

        Eigen::Vector3d rpy_init;
        Eigen::Vector3d rpy_final;

        state_vector_t x_init;
        state_vector_t x_final;

        double m_scale, r_scale;

        bool constrain_initial_final;

        void loadFromFile(const std::string &path);

        void nondimensionalize();

        void redimensionalize();
    } p;
};

} // namespace rocketHover