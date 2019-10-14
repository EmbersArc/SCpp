#include <cppad/cppad.hpp>

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
class SystemDynamics
{
public:
    using state_vector_t = Eigen::Matrix<double, STATE_DIM, 1>;
    using state_matrix_t = Eigen::Matrix<double, STATE_DIM, STATE_DIM>;
    using input_vector_t = Eigen::Matrix<double, INPUT_DIM, 1>;
    using input_matrix_t = Eigen::Matrix<double, INPUT_DIM, INPUT_DIM>;
    using control_matrix_t = Eigen::Matrix<double, STATE_DIM, INPUT_DIM>;
    using feedback_matrix_t = Eigen::Matrix<double, INPUT_DIM, STATE_DIM>;
    using param_vector_t = Eigen::Matrix<double, PARAM_DIM, 1>;
    using dynamic_vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using dynamic_matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using dynamic_vector_map_t = Eigen::Map<dynamic_vector_t>;

    using state_vector_v_t = std::vector<state_vector_t>;
    using input_vector_v_t = std::vector<input_vector_t>;
    using state_matrix_v_t = std::vector<state_matrix_t>;
    using control_matrix_v_t = std::vector<control_matrix_t>;

    using scalar_t = double;
    using scalar_ad_t = CppAD::AD<scalar_t>;

    using state_vector_ad_t = Eigen::Matrix<scalar_ad_t, STATE_DIM, 1>;
    using state_matrix_ad_t = Eigen::Matrix<scalar_ad_t, STATE_DIM, STATE_DIM>;
    using input_vector_ad_t = Eigen::Matrix<scalar_ad_t, INPUT_DIM, 1>;
    using control_matrix_ad_t = Eigen::Matrix<scalar_ad_t, STATE_DIM, INPUT_DIM>;
    using dynamic_vector_ad_t = Eigen::Matrix<scalar_ad_t, Eigen::Dynamic, 1>;
    using domain_vector_ad_t = Eigen::Matrix<scalar_ad_t, STATE_DIM + INPUT_DIM, 1>;
    using param_vector_ad_t = Eigen::Matrix<scalar_ad_t, PARAM_DIM, 1>;

    /**
     * @brief Initialize the model by compiling the dynamics functions
     * 
     */
    void initializeModel();

    /**
     * @brief Update model parameters.
     * 
     * @param param 
     */
    void updateModelParameters(param_vector_t param);

    /**
     * @brief The state derivative function. Has to be implemented by the derived class. All types have to be scalar_ad_t.
     * 
     * @param x 
     * @param u 
     * @param f 
     */
    virtual void systemFlowMap(
        const state_vector_ad_t &x,
        const input_vector_ad_t &u,
        const param_vector_ad_t &par,
        state_vector_ad_t &f) = 0;

    /**
     * @brief Compute the state derivative f(x,u)
     * 
     * @param x 
     * @param u 
     * @param f 
     */
    void computef(const state_vector_t &x, const input_vector_t &u, state_vector_t &f);

    /**
     * @brief Compute the state and control Jacobians A(x,u) and B(x,u)
     * 
     * @param x 
     * @param u 
     * @param A 
     * @param B 
     */
    void computeJacobians(const state_vector_t &x, const input_vector_t &u, state_matrix_t &A, control_matrix_t &B);

private:
    // CppAD function
    CppAD::ADFun<scalar_t> f_;
    bool initialized = false;
};

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
void SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::initializeModel()
{
    if (initialized)
    {
        return;
    }

    dynamic_vector_ad_t x(STATE_DIM + INPUT_DIM);
    dynamic_vector_ad_t param(PARAM_DIM);
    x.setOnes();
    param.setOnes();

    // start recording
    CppAD::Independent(x, 0, false, param);

    const state_vector_ad_t &state = x.head<STATE_DIM>();
    const input_vector_ad_t &input = x.tail<INPUT_DIM>();

    state_vector_ad_t dx;
    systemFlowMap(state, input, param, dx);

    // store operation sequence in x' = f(x) and stop recording
    f_ = CppAD::ADFun<scalar_t>(x, dynamic_vector_ad_t(dx));
    f_.optimize();

    initialized = true;
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
void SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::updateModelParameters(param_vector_t param)
{
    f_.new_dynamic(dynamic_vector_t(param));
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
void SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::computef(const state_vector_t &x,
                                                               const input_vector_t &u,
                                                               state_vector_t &f)
{
    assert(initialized);
    dynamic_vector_t input(STATE_DIM + INPUT_DIM, 1);
    input << x, u;
    dynamic_vector_map_t f_map(f.data(), STATE_DIM);

    f_map << f_.Forward(0, input);
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
void SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::computeJacobians(const state_vector_t &x,
                                                                       const input_vector_t &u,
                                                                       state_matrix_t &A,
                                                                       control_matrix_t &B)
{
    assert(initialized);
    dynamic_vector_t input(STATE_DIM + INPUT_DIM, 1);
    input << x, u;
    Eigen::Matrix<double, STATE_DIM, STATE_DIM + INPUT_DIM, Eigen::RowMajor> J;
    dynamic_vector_map_t J_map(J.data(), (STATE_DIM + INPUT_DIM) * STATE_DIM);

    J_map << f_.Jacobian(input);

    A = J.template leftCols<STATE_DIM>();
    B = J.template rightCols<INPUT_DIM>();
}
