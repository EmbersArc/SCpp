#pragma once

#include "activeModel.hpp"
#include "ecosWrapper.hpp"
#include "parameterServer.hpp"

namespace scpp
{

class SCvxRHAlgorithm
{
public:
    /**
     * @brief Construct a new SC solver.
     * 
     * @param model     The system model.
     */
    explicit SCvxRHAlgorithm(Model::ptr_t model);

    /**
     * @brief Initializes the algorithm. Has to be called before solving the problem.
     * 
     */
    void initialize();

    /**
     * @brief Solves the system.
     * 
     * @param warm_start    Whether to reuse the last computed trajectory.
     */
    void solve(bool warm_start = false);

    /**
     * @brief Get the solution variables object.
     * 
     */
    void getSolution(TrajectoryData &trajectory) const;

    /**
     * @brief Get the solution from each iteration
     * 
     */
    void getAllSolutions(std::vector<TrajectoryData> &all_trajectories) const;

    void setInitialState(const Model::state_vector_t &x);

    void setFinalState(const Model::state_vector_t &x);

private:
    /**
     * @brief Saves solution indices for performance.
     * 
     */
    void cacheIndices();

    /**
     * @brief Reads the solution variables X, U.
     * 
     */
    void readSolution();

    /**
     * @brief Loads the parameters from the configuration file.
     * 
     */
    void loadParameters();

    /**
     * @brief Performs a Successive Convexification iteration.
     * 
     */
    bool iterate();

    /**
     * @brief Get the nonlinear cost by integrating the dynamics
     * 
     * @return double 
     */
    double getNonlinearCost();

    void resetTrajectory();

    size_t K;

    Model::ptr_t model;

    bool interpolate_input;

    bool nondimensionalize;
    double alpha;
    double beta;
    double rho_0;
    double rho_1;
    double rho_2;
    double change_threshold;

    double trust_region;
    double initial_trust_region;
    double weight_virtual_control;
    size_t max_iterations;
    std::optional<double> last_nonlinear_cost;

    DiscretizationData dd;

    TrajectoryData td;

    Model::state_vector_t state_weights;
    Model::input_vector_t input_weights;

    std::vector<TrajectoryData> all_td;

    Eigen::MatrixXi X_indices;
    Eigen::MatrixXi U_indices;

    op::SecondOrderConeProgram socp;

    std::unique_ptr<EcosWrapper> solver;
};

} // namespace scpp