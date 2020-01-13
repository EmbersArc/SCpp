#pragma once

#include "activeModel.hpp"
#include "socpSolver.hpp"
#include "SCProblem.hpp"
#include "parameterServer.hpp"

namespace scpp
{

class SCAlgorithm
{
public:
    /**
     * @brief Construct a new SC solver.
     * 
     * @param model     The system model.
     */
    explicit SCAlgorithm(Model::ptr_t model);

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
     * @param X     The state trajectory.
     * @param U     The input trajectory.
     * @param t     The final time.
     */
    void getSolution(trajectory_data_t &trajectory) const;

    /**
     * @brief Get the solution from each iteration
     * 
     */
    void getAllSolutions(std::vector<trajectory_data_t> &all_trajectories);

private:
    /**
     * @brief Reads the solution variables X, U and sigma.
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
     * @return true     If converged.
     * @return false    If not converged.
     */
    bool iterate();

    size_t K;

    Model::ptr_t model;

    bool free_final_time;
    bool interpolate_input;

    bool nondimensionalize;
    double weight_time;
    double weight_trust_region_time;
    double weight_trust_region_trajectory;
    double weight_virtual_control;
    double trust_region_factor;
    double nu_tol;
    double delta_tol;
    size_t max_iterations;

    discretization_data_t dd;

    trajectory_data_t td;
    std::vector<trajectory_data_t> all_td;

    size_t sigma_index = 0;
    Eigen::MatrixXi X_indices;
    Eigen::MatrixXi U_indices;

    op::SecondOrderConeProgram socp;

    std::unique_ptr<op::Solver> solver;
};

} // namespace scpp