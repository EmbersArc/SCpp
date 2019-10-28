#pragma once

#define CODEGEN true

#if CODEGEN
#include <cppad/cg.hpp>
#else
#include <cppad/cppad.hpp>
#endif

#include "timing.hpp"

namespace scpp
{

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

#if CODEGEN
    using scalar_t = CppAD::cg::CG<double>;
#else
    using scalar_t = double;
#endif

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
    void computeJacobians(const state_vector_t &x,
                          const input_vector_t &u,
                          state_matrix_t &A,
                          control_matrix_t &B);

    struct DiscretizationData
    {
        state_matrix_v_t A;
        control_matrix_v_t B;
        control_matrix_v_t C;
        state_vector_v_t s;
        state_vector_v_t z;

        void initialize(size_t K, bool interpolate_input, bool free_final_time);

        bool interpolatedInput() const;

        bool variableTime() const;

        size_t n_X() const;
        size_t n_U() const;
    };

    struct TrajectoryData
    {
        state_vector_v_t X;
        input_vector_v_t U;
        double t;

        void initialize(size_t K, bool interpolate_input);

        size_t n_X() const;
        size_t n_U() const;
    };

private:
    CppAD::ADFun<scalar_t> f_;
#if CODEGEN
    std::unique_ptr<CppAD::cg::DynamicLib<double>> dynamicLib;
    std::unique_ptr<CppAD::cg::GenericModel<double>> model;
    param_vector_t current_parameters;
#endif

    bool initialized = false;
    bool parameters_set = false;
};

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
void SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::initializeModel()
{
    if (initialized)
    {
        return;
    }

#if CODEGEN
    dynamic_vector_ad_t x(STATE_DIM + INPUT_DIM + PARAM_DIM);
    x.setOnes();
    CppAD::Independent(x, 0, false);

    const state_vector_ad_t &state = x.segment<STATE_DIM>(0);
    const input_vector_ad_t &input = x.segment<INPUT_DIM>(STATE_DIM);
    const param_vector_ad_t &param = x.segment<PARAM_DIM>(STATE_DIM + INPUT_DIM);

    state_vector_ad_t dx;
    systemFlowMap(state, input, param, dx);

    f_ = CppAD::ADFun<scalar_t>(x, dynamic_vector_ad_t(dx));
    f_.optimize();

    CppAD::cg::ModelCSourceGen<double> cgen(f_, "model");
    cgen.setCreateForwardZero(true);
    cgen.setCreateJacobian(true);
    CppAD::cg::ModelLibraryCSourceGen<double> libcgen(cgen);

    // compile source code
    CppAD::cg::DynamicModelLibraryProcessor<double> p(libcgen);

    CppAD::cg::GccCompiler<double> compiler;
    compiler.addCompileFlag("-O3");
    dynamicLib = p.createDynamicLibrary(compiler);

    model = dynamicLib->model("model");
#else
    CppAD::thread_alloc::hold_memory(true);

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
#endif

    initialized = true;
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
void SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::updateModelParameters(param_vector_t param)
{
#if CODEGEN
    current_parameters = param;
#else
    f_.new_dynamic(dynamic_vector_t(param));
#endif
    parameters_set = true;
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
void SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::computef(const state_vector_t &x,
                                                               const input_vector_t &u,
                                                               state_vector_t &f)
{
    assert(initialized and parameters_set);

#if CODEGEN
    dynamic_vector_t input(STATE_DIM + INPUT_DIM + PARAM_DIM);
    input << x, u, current_parameters;

    CppAD::cg::ArrayView<const double> input_view(input.data(), input.size());
    CppAD::cg::ArrayView<double> f_view(f.data(), f.size());

    model->ForwardZero(input_view, f_view);
#else
    dynamic_vector_t input(STATE_DIM + INPUT_DIM);
    input << x, u;
    dynamic_vector_map_t f_map(f.data(), STATE_DIM);

    f_map << f_.Forward(0, input);
#endif
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
void SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::computeJacobians(const state_vector_t &x,
                                                                       const input_vector_t &u,
                                                                       state_matrix_t &A,
                                                                       control_matrix_t &B)
{
    assert(initialized and parameters_set);

#if CODEGEN
    dynamic_vector_t input(STATE_DIM + INPUT_DIM + PARAM_DIM);
    input << x, u, current_parameters;

    using full_jacobian_t = Eigen::Matrix<double, STATE_DIM, STATE_DIM + INPUT_DIM + PARAM_DIM, Eigen::RowMajor>;
    full_jacobian_t J;

    CppAD::cg::ArrayView<const double> input_view(input.data(), input.size());
    CppAD::cg::ArrayView<double> J_view(J.data(), J.size());

    model->Jacobian(input_view, J_view);
#else
    dynamic_vector_t input(STATE_DIM + INPUT_DIM);
    input << x, u;
    Eigen::Matrix<double, STATE_DIM, STATE_DIM + INPUT_DIM, Eigen::RowMajor> J;
    dynamic_vector_map_t J_map(J.data(), J.size());

    J_map << f_.Jacobian(input);
#endif

    A = J.template block<STATE_DIM, STATE_DIM>(0, 0);
    B = J.template block<STATE_DIM, INPUT_DIM>(0, STATE_DIM);
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
void SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::DiscretizationData::initialize(size_t K, bool interpolate_input, bool free_final_time)
{
    A.resize(K - 1);
    B.resize(K - 1);
    if (interpolate_input)
    {
        C.resize(K - 1);
    }
    if (free_final_time)
    {
        s.resize(K - 1);
    }
    z.resize(K - 1);
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
void SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::TrajectoryData::initialize(size_t K, bool interpolate_input)
{
    X.resize(K);
    U.resize(interpolate_input ? K : K - 1);
    t = 0.;
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
bool SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::DiscretizationData::interpolatedInput() const
{
    return not C.empty();
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
bool SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::DiscretizationData::variableTime() const
{
    return not s.empty();
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
size_t SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::DiscretizationData::n_X() const
{
    return A.size();
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
size_t SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::DiscretizationData::n_U() const
{
    return B.size();
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
size_t SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::TrajectoryData::n_X() const
{
    return X.size();
}

template <size_t STATE_DIM, size_t INPUT_DIM, size_t PARAM_DIM>
size_t SystemDynamics<STATE_DIM, INPUT_DIM, PARAM_DIM>::TrajectoryData::n_U() const
{
    return U.size();
}

} // namespace scpp