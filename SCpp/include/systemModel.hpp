#pragma once
#define JIT true

#include <Eigen/Dense>
#include <Eigen/StdVector>

#if JIT
#include <cppad/cg.hpp>
#else
#include <cppad/cppad.hpp>
#endif

#include "optimizationProblem.hpp"

template <size_t STATE_DIM, size_t INPUT_DIM>
class SystemModel
{
  public:
    enum : size_t
    {
        state_dim_ = STATE_DIM,
        input_dim_ = INPUT_DIM,
    };

    typedef Eigen::Matrix<double, STATE_DIM, 1> state_vector_t;
    typedef Eigen::Matrix<double, STATE_DIM, STATE_DIM> state_matrix_t;
    typedef Eigen::Matrix<double, INPUT_DIM, 1> input_vector_t;
    typedef Eigen::Matrix<double, STATE_DIM, INPUT_DIM> control_matrix_t;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dynamic_vector_t;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dynamic_matrix_t;
    typedef Eigen::Map<dynamic_vector_t> dynamic_vector_map_t;

#if JIT
    typedef CppAD::cg::CG<double> scalar_t;
#else
    typedef double scalar_t;
#endif

    typedef CppAD::AD<scalar_t> scalar_ad_t;

    typedef Eigen::Matrix<scalar_ad_t, STATE_DIM, 1> state_vector_ad_t;
    typedef Eigen::Matrix<scalar_ad_t, STATE_DIM, STATE_DIM> state_matrix_ad_t;
    typedef Eigen::Matrix<scalar_ad_t, INPUT_DIM, 1> input_vector_ad_t;
    typedef Eigen::Matrix<scalar_ad_t, STATE_DIM, INPUT_DIM> control_matrix_ad_t;
    typedef Eigen::Matrix<scalar_ad_t, Eigen::Dynamic, 1> dynamic_vector_ad_t;
    typedef Eigen::Matrix<scalar_ad_t, STATE_DIM + INPUT_DIM, 1> domain_vector_ad_t;

    typedef vector<state_vector_t, Eigen::aligned_allocator<state_vector_t>> state_vector_v_t;
    typedef vector<input_vector_t, Eigen::aligned_allocator<input_vector_t>> input_vector_v_t;
    typedef vector<state_matrix_t, Eigen::aligned_allocator<state_matrix_t>> state_matrix_v_t;
    typedef vector<control_matrix_t, Eigen::aligned_allocator<control_matrix_t>> control_matrix_v_t;

    /**
     * @brief Construct a new System Model object
     * 
     */
    SystemModel(){};

    /**
     * @brief Initialize the model by compiling the dynamics functions
     * 
     * @tparam STATE_DIM
     * @tparam INPUT_DIM 
     */
    void initializeModel();

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

    /**
     * @brief Function to initialize the trajectory of a derived model. Has to be implemented by the derived class.
     * 
     * @param X 
     * @param U 
     */
    virtual void initializeTrajectory(Eigen::MatrixXd &X,
                                      Eigen::MatrixXd &U) = 0;

    /**
     * @brief Function to add constraints of a derived model. Has to be implemented by the derived class.
     * 
     * @param socp 
     * @param X0 
     * @param U0 
     */
    virtual void addApplicationConstraints(
        op::SecondOrderConeProgram &socp,
        Eigen::MatrixXd &X0,
        Eigen::MatrixXd &U0) = 0;

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
        state_vector_ad_t &f) = 0;

    /**
     * @brief Function to remove mass and length dimensions from all function parameters.
     * 
     */
    virtual void nondimensionalize(){};
    /**
     * @brief Function to add mass and length dimensions to state and input trajectory.
     * 
     * @param X 
     * @param U 
     */
    virtual void redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                            Eigen::MatrixXd &U){};
    double r_scale = 1;
    double m_scale = 1;

    state_vector_t x_init;
    state_vector_t x_final;

    /**
     * @brief Get the final time guess
     * 
     * @return double The flight time guess
     */
    virtual double getFinalTimeGuess() = 0;

  private:
    /**
   * @brief Calculates the state derivative for AD. Uses the systemFlowMap() of the derived class.
   * 
   * @param tapedInput 
   * @param f 
   */
    void systemFlowMapAD(
        const dynamic_vector_ad_t &tapedInput,
        dynamic_vector_ad_t &f);

    // CppAD function
    CppAD::ADFun<scalar_t> f_;

#if JIT
    // Dynamic library and prepared model instances
    std::unique_ptr<CppAD::cg::DynamicLib<double>> dynamicLib_;
    std::unique_ptr<CppAD::cg::GenericModel<double>> model_;
#endif
};

template <size_t STATE_DIM, size_t INPUT_DIM>
void SystemModel<STATE_DIM, INPUT_DIM>::initializeModel()
{
    dynamic_vector_ad_t x(STATE_DIM + INPUT_DIM);
    x.setRandom();

    // start recording
    CppAD::Independent(x);

    dynamic_vector_ad_t dxFixed;
    systemFlowMapAD(x, dxFixed);

    // store operation sequence in x' = f(x) and stop recording
    f_ = CppAD::ADFun<scalar_t>(x, dxFixed);
    f_.optimize();

#if JIT
    // generate source code
    CppAD::cg::ModelCSourceGen<double> cgen(f_, "model");
    cgen.setCreateJacobian(true);
    cgen.setCreateHessian(true);
    cgen.setCreateForwardZero(true);
    CppAD::cg::ModelLibraryCSourceGen<double> libcgen(cgen);

    // compile source code
    CppAD::cg::GccCompiler<double> compiler;
    CppAD::cg::DynamicModelLibraryProcessor<double> processor(libcgen);
    compiler.setCompileFlags({"-O3", "-std=c11"});
    dynamicLib_ = processor.createDynamicLibrary(compiler);
    model_ = dynamicLib_->model("model");
#endif
}

template <size_t STATE_DIM, size_t INPUT_DIM>
void SystemModel<STATE_DIM, INPUT_DIM>::computef(const state_vector_t &x, const input_vector_t &u, state_vector_t &f)
{
    dynamic_vector_t input(STATE_DIM + INPUT_DIM, 1);
    input << x, u;
    dynamic_vector_map_t input_map(input.data(), STATE_DIM + INPUT_DIM);
    dynamic_vector_map_t f_map(f.data(), STATE_DIM);

#if JIT
    model_->ForwardZero(input_map, f_map);
#else
    f_map = f_.Forward(0, input);
#endif
}

template <size_t STATE_DIM, size_t INPUT_DIM>
void SystemModel<STATE_DIM, INPUT_DIM>::computeJacobians(const state_vector_t &x, const input_vector_t &u, state_matrix_t &A, control_matrix_t &B)
{
    dynamic_vector_t input(STATE_DIM + INPUT_DIM, 1);
    input << x, u;
    Eigen::Matrix<double, STATE_DIM, STATE_DIM + INPUT_DIM, Eigen::RowMajor> J;
    dynamic_vector_map_t input_map(const_cast<double *>(input.data()), STATE_DIM + INPUT_DIM);
    dynamic_vector_map_t J_map(J.data(), (STATE_DIM + INPUT_DIM) * STATE_DIM);

#if JIT
    model_->Jacobian(input_map, J_map);
#else
    J_map = f_.Jacobian(input);
#endif

    A = J.template leftCols<STATE_DIM>();
    B = J.template rightCols<INPUT_DIM>();
}

template <size_t STATE_DIM, size_t INPUT_DIM>
void SystemModel<STATE_DIM, INPUT_DIM>::systemFlowMapAD(
    const dynamic_vector_ad_t &tapedInput,
    dynamic_vector_ad_t &f)
{
    state_vector_ad_t x = tapedInput.template segment<STATE_DIM>(0);
    input_vector_ad_t u = tapedInput.template segment<INPUT_DIM>(STATE_DIM);

    state_vector_ad_t fFixed;
    systemFlowMap(x, u, fFixed);
    f = fFixed;
}
