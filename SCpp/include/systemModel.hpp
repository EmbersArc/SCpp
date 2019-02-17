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
     * @param X 
     * @param U 
     */
    virtual void initializeTrajectory(Eigen::MatrixXd &X,
                                      Eigen::MatrixXd &U) = 0;

    /**
     * @brief Function to add constraints of a derived model. Has to be implemented by the derived class.
     * 
     * @param socp 
     * @param X0 
     * @param U0 
     */
    virtual void addApplicationConstraints(
        op::SecondOrderConeProgram &socp,
        Eigen::MatrixXd &X0,
        Eigen::MatrixXd &U0) = 0;

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
     * @brief Get the final time guess
     * 
     * @return double The flight time guess
     */
    virtual double getFinalTimeGuess() = 0;

    /**
     * @brief Function to add mass and length dimensions to state and input trajectory.
     * 
     * @param X 
     * @param U 
     */
    virtual void redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                            Eigen::MatrixXd &U){};
    double r_scale = 1;
    double m_scale = 1;
};
