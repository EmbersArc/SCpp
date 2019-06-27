#include "activeModel.hpp"
#include "MPCProblem.hpp"
#include "discretization.hpp"
#include "ecosWrapper.hpp"
#include "parameterServer.hpp"

class MPCAlgorithm
{
public:
    /**
     * @brief Construct a new free Final Time Algorithm solver.
     * 
     * @param model     The system model.
     */
    explicit MPCAlgorithm(std::shared_ptr<Model> model);

    /**
     * @brief Initializes the algorithm. Has to be called before solving the problem.
     * 
     */
    void initialize();

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
    void loadParameters();

    ParameterServer param;

    size_t K;

    std::shared_ptr<Model> model;

    Model::state_matrix_t A;
    Model::control_matrix_t B;
    Model::state_vector_t z;

    Model::state_vector_t x_init;
    Model::state_vector_t x_final;

    Eigen::MatrixXd X;
    Eigen::MatrixXd U;

    Eigen::MatrixXi X_indices;
    Eigen::MatrixXi U_indices;

    op::SecondOrderConeProgram socp;

    std::unique_ptr<EcosWrapper> solver;
};