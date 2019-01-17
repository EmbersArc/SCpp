#pragma once

// #include <cppad/cg.hpp>
#include <cppad/cppad.hpp>

#include "optimizationProblem.hpp"

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM, size_t K>
class SystemModel
{
  public:
    typedef Eigen::Matrix<double, STATE_DIM, 1> state_vector_t;
    typedef Eigen::Matrix<double, STATE_DIM, STATE_DIM> state_matrix_t;
    typedef Eigen::Matrix<double, INPUT_DIM, 1> input_vector_t;
    typedef Eigen::Matrix<double, INPUT_DIM, INPUT_DIM> input_matrix_t;
    typedef Eigen::Matrix<double, STATE_DIM, INPUT_DIM> control_matrix_t;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dynamic_vector_t;

    typedef CppAD::AD<double> scalar_ad_t;
    // typedef CppAD::cg::CG<double> scalar_cg_ad_t;
    // typedef std::conditional<false, CppAD::AD<scalar_cg_ad_t>, CppAD::AD<default_scalar_ad_t>>::type scalar_ad_t;
    typedef Eigen::Matrix<scalar_ad_t, STATE_DIM, 1> state_vector_ad_t;
    typedef Eigen::Matrix<scalar_ad_t, STATE_DIM, STATE_DIM> state_matrix_ad_t;
    typedef Eigen::Matrix<scalar_ad_t, INPUT_DIM, 1> input_vector_ad_t;
    typedef Eigen::Matrix<scalar_ad_t, INPUT_DIM, INPUT_DIM> input_matrix_ad_t;
    typedef Eigen::Matrix<scalar_ad_t, STATE_DIM, INPUT_DIM> control_matrix_ad_t;
    typedef Eigen::Matrix<scalar_ad_t, Eigen::Dynamic, 1> dynamic_vector_ad_t;

    enum : size_t
    {
        state_dim_ = STATE_DIM,
        input_dim_ = INPUT_DIM,
    };

    SystemModel(){};
    void initializeModel();
    void computef(const state_vector_t &x, const input_vector_t &u, state_vector_t &f);
    void computeA(const state_vector_t &x, const input_vector_t &u, state_matrix_t &A);
    void computeB(const state_vector_t &x, const input_vector_t &u, control_matrix_t &B);

    template <typename T>
    void systemFlowMap(
        const Eigen::Matrix<T, STATE_DIM, 1> &x,
        const Eigen::Matrix<T, INPUT_DIM, 1> &u,
        Eigen::Matrix<T, STATE_DIM, 1> &f);

    void systemFlowMapAD(
        const dynamic_vector_ad_t &tapedInput,
        dynamic_vector_ad_t &f);

    void initializeTrajectory(Eigen::Matrix<double, STATE_DIM, K> &X,
                              Eigen::Matrix<double, INPUT_DIM, K> &U);

    virtual void addApplicationConstraints(
        optimization_problem::SecondOrderConeProgram &socp,
        Eigen::Matrix<double, STATE_DIM, K> &X0,
        Eigen::Matrix<double, INPUT_DIM, K> &U0) = 0;

  private:
    CppAD::ADFun<double> f_;
    vector<bool> sparsityA_, sparsityB_;
    vector<size_t> rowsA_, colsA_;
    vector<size_t> rowsB_, colsB_;
    CppAD::sparse_jacobian_work workA_, workB_;
};

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM, size_t K>
void SystemModel<Derived, STATE_DIM, INPUT_DIM, K>::initializeModel()
{
    // set up spasity
    {
        rowsA_.resize(STATE_DIM * STATE_DIM);
        colsA_.resize(STATE_DIM * STATE_DIM);
        rowsB_.resize(STATE_DIM * INPUT_DIM);
        colsB_.resize(STATE_DIM * INPUT_DIM);
        size_t countA = 0;
        size_t countB = 0;
        for (size_t i = 0; i < STATE_DIM + INPUT_DIM; i++)
        {
            for (size_t j = 0; j < STATE_DIM; j++)
            {
                if (i < STATE_DIM)
                {
                    rowsA_[countA] = j;
                    colsA_[countA] = i;
                    countA++;
                }
                else
                {
                    rowsB_[countB] = j;
                    colsB_[countB] = i;
                    countB++;
                }
            }
        }
    }

    sparsityA_.resize((STATE_DIM + INPUT_DIM) * STATE_DIM);
    sparsityB_.resize((STATE_DIM + INPUT_DIM) * STATE_DIM);
    for (size_t i = 0; i < sparsityA_.size(); i++)
    {
        if (i % (STATE_DIM + INPUT_DIM) < STATE_DIM)
        {
            sparsityA_[i] = true;
        }
        else
        {
            sparsityB_[i] = true;
        }
    }

    // record the model
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
        CppAD::ADFun<double> f(x, dxFixed);

        f.optimize();

        f_ = f;
    }
}

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM, size_t K>
void SystemModel<Derived, STATE_DIM, INPUT_DIM, K>::computef(const state_vector_t &x, const input_vector_t &u, state_vector_t &f)
{
    dynamic_vector_t input(STATE_DIM + INPUT_DIM);
    input << x, u;

    f = f_.Forward(0, input);
}

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM, size_t K>
void SystemModel<Derived, STATE_DIM, INPUT_DIM, K>::computeA(const state_vector_t &x, const input_vector_t &u, state_matrix_t &A)
{
    dynamic_vector_t input(STATE_DIM + INPUT_DIM);
    input << x, u;

    dynamic_vector_t jac(STATE_DIM * STATE_DIM);
    f_.SparseJacobianForward(input, sparsityA_, rowsA_, colsA_, jac, workA_);

    Eigen::Map<state_matrix_t> out(jac.data());

    A = out;
}

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM, size_t K>
void SystemModel<Derived, STATE_DIM, INPUT_DIM, K>::computeB(const state_vector_t &x, const input_vector_t &u, control_matrix_t &B)
{
    dynamic_vector_t input(STATE_DIM + INPUT_DIM);
    input << x, u;

    dynamic_vector_t jac(STATE_DIM * INPUT_DIM);
    f_.SparseJacobianForward(input, sparsityB_, rowsB_, colsB_, jac, workB_);

    Eigen::Map<control_matrix_t> out(jac.data());

    B = out;
}

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM, size_t K>
template <typename T>
void SystemModel<Derived, STATE_DIM, INPUT_DIM, K>::systemFlowMap(
    const Eigen::Matrix<T, STATE_DIM, 1> &x,
    const Eigen::Matrix<T, INPUT_DIM, 1> &u,
    Eigen::Matrix<T, STATE_DIM, 1> &f)
{
    throw std::runtime_error("systemFlowMap() method should be implemented by the derived class.");
}

template <class Derived, size_t STATE_DIM, size_t INPUT_DIM, size_t K>
void SystemModel<Derived, STATE_DIM, INPUT_DIM, K>::systemFlowMapAD(
    const dynamic_vector_ad_t &tapedInput,
    dynamic_vector_ad_t &f)
{
    state_vector_ad_t x = tapedInput.segment(0, STATE_DIM);
    input_vector_ad_t u = tapedInput.segment(STATE_DIM, INPUT_DIM);

    state_vector_ad_t fFixed;
    static_cast<Derived *>(this)->template systemFlowMap<scalar_ad_t>(x, u, fFixed);
    f = fFixed;
}