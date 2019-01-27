#pragma once

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <cppad/cg.hpp>

#include "optimizationProblem.hpp"

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM>
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

    typedef CppAD::cg::CG<double> scalar_cg_ad_t;
    typedef CppAD::AD<scalar_cg_ad_t> scalar_ad_t;
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

    // Default constructor
    SystemModel(){};

    // Set up and compile dynamics function
    void initializeModel();

    // Computes the state derivative
    void computef(const state_vector_t &x, const input_vector_t &u, state_vector_t &f, size_t modelNum = 0);
    // Computes state and control Jacobians
    void computeJacobians(const state_vector_t &x, const input_vector_t &u, state_matrix_t &A, control_matrix_t &B = 0);

    // Function to initialize the trajectory of a derived model. Has to be implemented by the derived class.
    virtual void initializeTrajectory(Eigen::MatrixXd &X,
                                      Eigen::MatrixXd &U) = 0;

    // Function to add constraints of a derived model. Has to be implemented by the derived class.
    virtual void addApplicationConstraints(
        op::SecondOrderConeProgram &socp,
        Eigen::MatrixXd &X0,
        Eigen::MatrixXd &U0) = 0;

    // The state derivative function. Has to be implemented by the derived class.
    virtual void systemFlowMap(
        const state_vector_ad_t &x,
        const input_vector_ad_t &u,
        state_vector_ad_t &f) = 0;

    virtual void getStateWeightVector(state_vector_t &w) = 0;
    virtual void getInputWeightVector(input_vector_t &w) = 0;
    virtual void nondimensionalize() = 0;

  private:
    // Calculates the state derivative for AD. Uses the systemFlowMap() of the derived class.
    void systemFlowMapAD(
        const dynamic_vector_ad_t &tapedInput,
        dynamic_vector_ad_t &f);

    // Dynamic library and prepared model instances
    std::unique_ptr<CppAD::cg::DynamicLib<double>> dynamicLib_;
    std::unique_ptr<CppAD::cg::GenericModel<double>> model_;
};

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM>
void SystemModel<Derived, STATE_DIM, INPUT_DIM>::initializeModel()
{
    // input vector
    dynamic_vector_ad_t x(STATE_DIM + INPUT_DIM);
    x.setRandom();

    // declare x as independent
    CppAD::Independent(x);

    // create fixed size types since CT uses fixed size types
    dynamic_vector_ad_t dxFixed;

    systemFlowMapAD(x, dxFixed);

    // store operation sequence in f: x -> dx and stop recording
    CppAD::ADFun<scalar_cg_ad_t> f(x, dxFixed);

    f.optimize();

    // generates source code
    CppAD::cg::ModelCSourceGen<double> cgen(f, "model");
    cgen.setCreateJacobian(true);
    cgen.setCreateForwardZero(true);
    CppAD::cg::ModelLibraryCSourceGen<double> libcgen(cgen);

    // compile source code
    CppAD::cg::GccCompiler<double> compiler;
    CppAD::cg::DynamicModelLibraryProcessor<double> processor(libcgen);
    compiler.setCompileFlags({"-O3", "-std=c11"});
    dynamicLib_ = processor.createDynamicLibrary(compiler);

    // save to files
    // CppAD::cg::SaveFilesModelLibraryProcessor<double> processorSave(libcgen);
    // processorSave.saveSources();

    model_ = dynamicLib_->model("model");
}

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM>
void SystemModel<Derived, STATE_DIM, INPUT_DIM>::computef(const state_vector_t &x, const input_vector_t &u, state_vector_t &f, size_t modelNum)
{
    dynamic_vector_t input(STATE_DIM + INPUT_DIM, 1);
    input << x, u;
    dynamic_vector_map_t input_map(input.data(), STATE_DIM + INPUT_DIM);
    dynamic_vector_map_t f_map(f.data(), STATE_DIM);

    model_->ForwardZero(input_map, f_map);
}

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM>
void SystemModel<Derived, STATE_DIM, INPUT_DIM>::computeJacobians(const state_vector_t &x, const input_vector_t &u, state_matrix_t &A, control_matrix_t &B)
{
    dynamic_vector_t input(STATE_DIM + INPUT_DIM, 1);
    input << x, u;
    dynamic_vector_map_t input_map(const_cast<double *>(input.data()), STATE_DIM + INPUT_DIM);
    Eigen::Matrix<double, STATE_DIM, STATE_DIM + INPUT_DIM, Eigen::RowMajor> J;
    dynamic_vector_map_t J_map(J.data(), (STATE_DIM + INPUT_DIM) * STATE_DIM);

    model_->Jacobian(input_map, J_map);
    A = J.template leftCols<STATE_DIM>();
    B = J.template rightCols<INPUT_DIM>();
}

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM>
void SystemModel<Derived, STATE_DIM, INPUT_DIM>::systemFlowMapAD(
    const dynamic_vector_ad_t &tapedInput,
    dynamic_vector_ad_t &f)
{
    state_vector_ad_t x = tapedInput.template segment<STATE_DIM>(0);
    input_vector_ad_t u = tapedInput.template segment<INPUT_DIM>(STATE_DIM);

    state_vector_ad_t fFixed;
    systemFlowMap(x, u, fFixed);
    f = fFixed;
}
