#pragma once

#include <Eigen/Dense>

#include "epigraph.hpp"
#include "systemDynamics.hpp"
#include "trajectoryData.hpp"
#include "discretizationData.hpp"

namespace scpp
{

template <typename DERIVED, size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
class SystemModel : public SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>
{
public:
    using BASE = SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>;

    using trajectory_data_t = TrajectoryData<BASE>;
    using discretization_data_t = DiscretizationData<BASE>;

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

    virtual ~SystemModel(){};

    /**
     * @brief Function to initialize the trajectory of a derived model. Has to be implemented by the derived class. Only required for SC models,
     * 
     */
    virtual void getInitializedTrajectory(trajectory_data_t &)
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
    virtual void addApplicationConstraints(std::shared_ptr<cvx::OptimizationProblem>,
                                           state_vector_v_t &,
                                           input_vector_v_t &){};

    /**
     * @brief Gets the new parameters in the system flow map.
     * 
     */
    virtual void getNewModelParameters(param_vector_t &){};

    /**
     * @brief Updates the parameters in the system flow map.
     * 
     */
    void updateModelParameters()
    {
        param_vector_t model_params;
        getNewModelParameters(model_params);
        BASE::updateModelParameters(model_params);
    };

    /**
     * @brief Function to remove mass and length dimensions from all function parameters.
     * 
     */
    virtual void nondimensionalize()
    {
        throw std::runtime_error("nondimensionalize: This function has to be implemented by the derived class.");
    };

    /**
     * @brief Function to add mass and length dimensions to all function parameters.
     * 
     */
    virtual void redimensionalize()
    {
        throw std::runtime_error("redimensionalize: This function has to be implemented by the derived class.");
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
     */
    virtual void nondimensionalizeTrajectory(trajectory_data_t &)
    {
        throw std::runtime_error("nondimensionalizeTrajectory: This function has to be implemented by the derived class.");
    };

    /**
     * @brief Function to add mass and length dimensions to state and input trajectory.
     * 
     */
    virtual void redimensionalizeTrajectory(trajectory_data_t &)
    {
        throw std::runtime_error("redimensionalizeTrajectory: This function has to be implemented by the derived class.");
    };

    static const std::string getModelName()
    {
        return DERIVED::modelName;
    };

    void setParameterFolder(const std::string &path)
    {
        param_folder_path = path;
    };

    const std::string getParameterFolder() const
    {
        return param_folder_path + getModelName();
    };

private:
    std::string param_folder_path = "../scpp_models/config/";
};

} // namespace scpp