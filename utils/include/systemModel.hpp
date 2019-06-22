#pragma once

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "optimizationProblem.hpp"
#include "systemDynamics.hpp"

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
class SystemModel : public SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>
{
public:
    typedef SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM> BASE;

    using typename BASE::control_matrix_t;
    using typename BASE::control_matrix_v_t;
    using typename BASE::dynamic_matrix_t;
    using typename BASE::dynamic_vector_map_t;
    using typename BASE::dynamic_vector_t;
    using typename BASE::input_vector_t;
    using typename BASE::input_vector_v_t;
    using typename BASE::param_vector_t;
    using typename BASE::state_matrix_t;
    using typename BASE::state_matrix_v_t;
    using typename BASE::state_vector_t;
    using typename BASE::state_vector_v_t;

    using typename BASE::control_matrix_ad_t;
    using typename BASE::domain_vector_ad_t;
    using typename BASE::dynamic_vector_ad_t;
    using typename BASE::input_vector_ad_t;
    using typename BASE::param_vector_ad_t;
    using typename BASE::state_matrix_ad_t;
    using typename BASE::state_vector_ad_t;

    enum : size_t
    {
        state_dim = STATE_DIM,
        input_dim = INPUT_DIM,
        param_dim = PARAM_DIM,
    };

    /**
     * @brief Construct a new System Model object
     * 
     */
    SystemModel(){};

    /**
     * @brief Function to initialize the trajectory of a derived model. Has to be implemented by the derived class.
     * 
     * @param X     State trajectory.
     * @param U     Input trajectory.
     */
    virtual void getInitializedTrajectory(Eigen::MatrixXd &X,
                                          Eigen::MatrixXd &U)
    {
        throw std::runtime_error("getInitializedTrajectory: This function has to be implemented by the derived class.");
    };

    /**
     * @brief Function to add constraints of a model. Has to be implemented by the derived class.
     * 
     * @param socp  The SOCP.
     * @param X     Last state trajectory.
     * @param U     Last input trajectory.
     */
    virtual void addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                           Eigen::MatrixXd &X0,
                                           Eigen::MatrixXd &U0){};

    /**
     * @brief Function to remove mass and length dimensions from all function parameters.
     * 
     */
    virtual void nondimensionalize(){};

    /**
     * @brief Function to add mass and length dimensions to all function parameters.
     * 
     */
    virtual void redimensionalize(){};

    /**
     * @brief Get the time horizon or initial time guess.
     * 
     * @return double The time horizon or initial time guess
     */
    virtual void getTimeHorizon(double &sigma)
    {
        throw std::runtime_error("getTimeHorizon: This function has to be implemented by the derived class.");
    };

    virtual void getStateWeights(state_vector_t &intermediate, state_vector_t &terminal)
    {
        throw std::runtime_error("getStateWeights: This function has to be implemented by the derived class.");
    };

    virtual void getInputWeights(input_vector_t &intermediate)
    {
        throw std::runtime_error("getInputWeights: This function has to be implemented by the derived class.");
    };

    /**
     * @brief Function to remove mass and length dimensions from state and input trajectory.
     * 
     * @param X     State trajectory.
     * @param U     Input trajectory.
     */
    virtual void nondimensionalizeTrajectory(Eigen::MatrixXd &X,
                                             Eigen::MatrixXd &U){};

    /**
     * @brief Function to add mass and length dimensions to state and input trajectory.
     * 
     * @param X     State trajectory.
     * @param U     Input trajectory.
     */
    virtual void redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                            Eigen::MatrixXd &U){};
};
