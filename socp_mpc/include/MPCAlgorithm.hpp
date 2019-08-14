#include "MPCProblem.hpp"
#include "activeModel.hpp"
#include "discretization.hpp"
#include "ecosWrapper.hpp"
#include "parameterServer.hpp"

namespace mpc
{

class MPCAlgorithm
{
public:
    /**
     * @brief Construct a new free Final Time Algorithm solver.
     *
     * @param model     The system model.
     */
    explicit MPCAlgorithm(std::shared_ptr<Model> model, const std::string parameter_path = "");

    /**
     * @brief Initializes the algorithm. Has to be called before solving the problem.
     *
     */
    void initialize(bool constant_dynamics = false,
                    bool intermediate_cost_active = true);

    /**
     * @brief Discretizes the system.
     *
     */
    void discretize(const Model::state_vector_t &x, const Model::input_vector_t &u);

    /**
     * @brief Returns the number of discretization steps.
     *
     */
    void getTimeSteps(size_t &K);

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
    void getSolution(Model::dynamic_matrix_t &X, Model::dynamic_matrix_t &U);

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
    void loadParameters(const std::string &path);

    size_t K;
    double dt;
    bool nondimensionalize;

    std::shared_ptr<Model> model;

    Model::state_matrix_t A;
    Model::control_matrix_t B;
    Model::state_vector_t z;

    Model::state_vector_t state_weights_intermediate;
    Model::state_vector_t state_weights_terminal;
    Model::input_vector_t input_weights;

    Model::state_vector_t x_init;
    Model::state_vector_t x_final;

    Eigen::MatrixXd X;
    Eigen::MatrixXd U;

    Eigen::MatrixXi X_indices;
    Eigen::MatrixXi U_indices;

    op::SecondOrderConeProgram socp;

    std::unique_ptr<EcosWrapper> solver;
};

} // namespace mpc