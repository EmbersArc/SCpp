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
class RocketHover : public SystemModel<STATE_DIM_, INPUT_DIM_, PARAM_DIM_>
{
public:
    RocketHover();

    static std::string getModelName()
    {
        return "RocketHover";
    }

    void systemFlowMap(
        const state_vector_ad_t &x,
        const input_vector_ad_t &u,
        const param_vector_ad_t &par,
        state_vector_ad_t &f) override;

    void getOperatingPoint(state_vector_t &x, input_vector_t &u);

    void addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                   Eigen::MatrixXd &X0,
                                   Eigen::MatrixXd &U0) override;

    void nondimensionalize() override;

    void redimensionalize() override;

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

        state_vector_t x_init;
        state_vector_t x_final;

        double m_scale, r_scale;

        void loadFromFile();

        void nondimensionalize();

        void redimensionalize();
    } p;
};

} // namespace rocketHover