#pragma once

#include <string>
#include <random>

#include "systemModel.hpp"
#include "ecosWrapper.hpp"
#include "parameterServer.hpp"

#include "rocketLanding3dDefinitions.hpp"

namespace rocketlanding3d
{

/**
 * @brief A 3D rocket landing model.
 * 
 */
class RocketLanding3D : public SystemModel<RocketLanding3D, STATE_DIM_, INPUT_DIM_, PARAM_DIM_>
{
public:
    RocketLanding3D();

    static std::string getModelName()
    {
        return "RocketLanding3D";
    }

    void systemFlowMap(
        const state_vector_ad_t &x,
        const input_vector_ad_t &u,
        const param_vector_ad_t &par,
        state_vector_ad_t &f) override;

    void getInitializedTrajectory(state_vector_v_t &X,
                                  input_vector_v_t &U) override;

    void addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                   state_vector_v_t &X0,
                                   input_vector_v_t &U0) override;

    void nondimensionalize() override;

    void redimensionalize() override;

    void nondimensionalizeTrajectory(state_vector_v_t &X,
                                     input_vector_v_t &U) override;

    void redimensionalizeTrajectory(state_vector_v_t &X,
                                    input_vector_v_t &U) override;

    void getNewModelParameters(param_vector_t &param) override;

    struct Parameters
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

        void randomizeInitialState();

        void loadFromFile();

        void nondimensionalize();

        void redimensionalize();

        void nondimensionalizeTrajectory(state_vector_v_t &X, input_vector_v_t &U) const;

        void redimensionalizeTrajectory(state_vector_v_t &X, input_vector_v_t &U) const;
    } p;
};

} // namespace rocketlanding3d