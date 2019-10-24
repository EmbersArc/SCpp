#pragma once

#include "activeModel.hpp"
#include "ecosWrapper.hpp"
#include "discretization.hpp"
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
    void getSolution(Model::state_vector_v_t &X, Model::input_vector_v_t &U, double &t);

    /**
     * @brief Get the solution from each iteration
     * 
     * @param X 
     * @param U 
     * @param t 
     */
    void getAllSolutions(std::vector<Model::state_vector_v_t> &X,
                         std::vector<Model::input_vector_v_t> &U,
                         std::vector<double> &t);

private:
    /**
     * @brief Saves solution indices for performance.
     * 
     */
    void cacheIndices();

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

    discretization::Data dd;

    double sigma;
    Model::state_vector_v_t X;
    Model::input_vector_v_t U;

    std::vector<Model::state_vector_v_t> all_X;
    std::vector<Model::input_vector_v_t> all_U;
    std::vector<double> all_times;

    size_t sigma_index = 0;
    Eigen::MatrixXi X_indices;
    Eigen::MatrixXi U_indices;

    op::SecondOrderConeProgram socp;

    std::unique_ptr<EcosWrapper> solver;
};

} // namespace scpp