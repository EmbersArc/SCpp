#pragma once

#include <string>
#include <random>

#include "systemModel.hpp"
#include "ecosWrapper.hpp"
#include "parameterServer.hpp"

#include "templateModelDefinitions.hpp"

namespace templateModel
{

/**
 * @brief A model template.
 * 
 */
class TemplateModel : public SystemModel<STATE_DIM_, INPUT_DIM_, PARAM_DIM_>
{
  public:
    TemplateModel();

    static std::string getModelName()
    {
        return "TemplateModel";
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

    void addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                   Eigen::MatrixXd &X0,
                                   Eigen::MatrixXd &U0) override;

    void nondimensionalize() override;

    void redimensionalize() override;

    void nondimensionalizeTrajectory(Eigen::MatrixXd &X,
                                     Eigen::MatrixXd &U) override;

    void redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                    Eigen::MatrixXd &U) override;

    void getNewModelParameters(param_vector_t &param);

    struct Parameters
    {
        state_vector_t x_init;
        state_vector_t x_final;
        double final_time_guess;

        double m_scale, r_scale;

        double m;

        void loadFromFile();

        void nondimensionalize();

        void redimensionalize();

        void nondimensionalizeTrajectory(Eigen::MatrixXd &X, Eigen::MatrixXd &U) const;

        void redimensionalizeTrajectory(Eigen::MatrixXd &X, Eigen::MatrixXd &U) const;
    } p;
};

} // namespace templateModel