#include "LQR.hpp"
#include "parameterServer.hpp"

namespace scpp
{

class LQRAlgorithm
{
public:
    /**
     * @brief Construct a new LQR solver.
     *
     * @param model     The system model.
     */
    explicit LQRAlgorithm(Model::ptr_t model);

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
     * @param weights 
     */
    void setStateWeights(const Model::state_vector_t &weights);

    /**
     * @brief Set the input weights
     * 
     * @param weights 
     */
    void setInputWeights(const Model::input_vector_t &weights);

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
    void getSolution(Model::input_vector_t &u);

private:
    /**
     * @brief Loads the parameters from the configuration file.
     *
     */
    void loadParameters();

    Model::feedback_matrix_t K;

    Model::ptr_t model;

    Model::state_vector_t x_eq;
    Model::input_vector_t u_eq;

    Model::state_matrix_t A;
    Model::control_matrix_t B;

    Model::state_matrix_t Q;
    Model::input_matrix_t R;

    std::optional<Model::input_vector_t> u;

    Model::state_vector_t state_weights;
    Model::input_vector_t input_weights;

    Model::state_vector_t x_init;
    Model::state_vector_t x_final;

    bool state_weights_set = false;
    bool input_weights_set = false;
    bool initialized = false;
};

} // namespace scpp