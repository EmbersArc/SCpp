#pragma once

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "fmt/format.h"

#include "optimizationProblem.hpp"
#include "systemDynamics.hpp"

template <typename DERIVED, size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
class SystemModel : public SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>
{
public:
    typedef SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM> BASE;

    using ptr_t = std::shared_ptr<DERIVED>;

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
     * @brief Function to initialize the trajectory of a derived model. Has to be implemented by the derived class. Only required for SC models,
     * 
     * @param X     State trajectory.
     * @param U     Input trajectory.
     */
    virtual void getInitializedTrajectory(Eigen::MatrixXd &,
                                          Eigen::MatrixXd &)
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
    virtual void addApplicationConstraints(op::SecondOrderConeProgram &,
                                           Eigen::MatrixXd &,
                                           Eigen::MatrixXd &){};

    /**
     * @brief Updates the parameters in the system flow map.
     * 
     */
    virtual void getNewModelParameters(param_vector_t &){};

    /**
     * @brief Function to remove mass and length dimensions from all function parameters.
     * 
     */
    virtual void nondimensionalize()
    {
        fmt::print("nondimensionalize: Function called without implementation.");
    };

    /**
     * @brief Function to add mass and length dimensions to all function parameters.
     * 
     */
    virtual void redimensionalize()
    {
        fmt::print("redimensionalize: Function called without implementation.");
    };

    /**
     * @brief Get the operating point of the system. Usually an equilibrium point for linearization.
     * 
     */
    virtual void getOperatingPoint(state_vector_t &, input_vector_t &)
    {
        throw std::runtime_error("getOperatingPoint: This function has to be implemented by the derived class.");
    };

    /**
     * @brief Function to remove mass and length dimensions from state and input trajectory.
     * 
     * @param X     State trajectory.
     * @param U     Input trajectory.
     */
    virtual void nondimensionalizeTrajectory(Eigen::MatrixXd &,
                                             Eigen::MatrixXd &)
    {
        fmt::print("nondimensionalizeTrajectory: Function called without implementation.");
    };

    /**
     * @brief Function to add mass and length dimensions to state and input trajectory.
     * 
     * @param X     State trajectory.
     * @param U     Input trajectory.
     */
    virtual void redimensionalizeTrajectory(Eigen::MatrixXd &,
                                            Eigen::MatrixXd &)
    {
        fmt::print("redimensionalizeTrajectory: Function called without implementation.");
    };

    void setParameterFolder(const std::string &path)
    {
        param_folder_path = path;
    }

    const std::string getParameterFolder()
    {
        return param_folder_path;
    }

private:
    std::string param_folder_path = "../socp_mpc/config";
};
