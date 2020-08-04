#pragma once

#include "activeModel.hpp"
#include "socpInterface.hpp"
#include "parameterServer.hpp"

namespace scpp
{

class MPCAlgorithm
{
public:
    /**
     * @brief Construct a new MPC solver.
     *
     * @param model     The system model.
     */
    explicit MPCAlgorithm(Model::ptr_t model);

    /**
     * @brief Initializes the algorithm. Has to be called before solving the problem.
     *
     */
    void initialize();

    /**
     * @brief  Sets a new initial state.
     *
     */
    void setInitialState(const Model::state_vector_t &x);

    /**
     * @brief  Sets a new desired state to track.
     *
     */
    void setFinalState(const Model::state_vector_t &x);

    /**
     * @brief Set the state weights
     * 
     * @param intermediate 
     * @param terminal 
     */
    void setStateWeights(const Model::state_vector_t &intermediate, const Model::state_vector_t &terminal);

    /**
     * @brief Set the input weights
     * 
     * @param intermediate 
     */
    void setInputWeights(const Model::input_vector_t &intermediate);

    /**
     * @brief Solves the system.
     *
     */
    void solve();

    /**
     * @brief Get the solution variables object.
     *
     * @param X     The state trajectory.
     * @param U     The input trajectory.
     */
    void getSolution(Model::state_vector_v_t &X, Model::input_vector_v_t &U);

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

    double time_horizon;
    size_t K;
    bool nondimensionalize;
    bool constant_dynamics;
    bool intermediate_cost_active;
    Model::state_vector_t state_weights_intermediate;
    Model::state_vector_t state_weights_terminal;
    Model::input_vector_t input_weights;

    Model::ptr_t model;

    Model::state_matrix_t A;
    Model::control_matrix_t B;
    Model::state_vector_t z;

    Model::state_vector_t x_init;
    Model::state_vector_t x_final;

    Model::state_vector_v_t X;
    Model::input_vector_v_t U;

    Eigen::MatrixXi X_indices;
    Eigen::MatrixXi U_indices;

    cvx::OptimizationProblem socp;

    std::unique_ptr<op::Solver> solver;

    bool state_weights_set = false;
    bool input_weights_set = false;
    bool initialized = false;
};

} // namespace scpp