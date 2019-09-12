#include "LQRAlgorithm.hpp"

namespace lqr
{

LQRAlgorithm::LQRAlgorithm(std::shared_ptr<Model> model, const std::string parameter_path) : model(model)
{
}

void LQRAlgorithm::initialize()
{
    assert(state_weights_set and input_weights_set);

    model->getOperatingPoint(x_eq, u_eq);
    model->computeJacobians(x_eq, u_eq, A, B);
    ComputeLQR(Q, R, A, B, K);

    initialized = true;
}

void LQRAlgorithm::solve()
{
    assert(initialized);
    const Model::state_vector_t state_error = x_init - x_final;
    u = -K * state_error + u_eq;
}

void LQRAlgorithm::setInitialState(const Model::state_vector_t &x)
{
    x_init = x;
}

void LQRAlgorithm::setFinalState(const Model::state_vector_t &x)
{
    x_final = x;
}

void LQRAlgorithm::setStateWeights(const Model::state_vector_t &weights)
{
    Q.diagonal() = weights;
    state_weights_set = true;
}

void LQRAlgorithm::setInputWeights(const Model::input_vector_t &weights)
{
    R.diagonal() = weights;
    input_weights_set = true;
}

void LQRAlgorithm::getSolution(Model::input_vector_t &u)
{
    u = this->u;
}

} // namespace lqr