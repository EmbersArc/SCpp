#include "LQRAlgorithm.hpp"

namespace scpp
{

LQRAlgorithm::LQRAlgorithm(Model::ptr_t model) : model(model)
{
    loadParameters();
}

void LQRAlgorithm::initialize()
{
    // print("[LQR] Starting controller for model '{}'.\n", Model::getModelName());

    assert(state_weights_set and input_weights_set);

    model->updateModelParameters();

    model->getOperatingPoint(x_eq, u_eq);
    model->computeJacobians(x_eq, u_eq, A, B);
    ComputeLQR(Q, R, A, B, K);

    initialized = true;
    // print("[LQR] Controller started.\n");
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
    Q.setZero();
    Q.diagonal() = weights;
    state_weights_set = true;
}

void LQRAlgorithm::setInputWeights(const Model::input_vector_t &weights)
{
    R.setZero();
    R.diagonal() = weights;
    input_weights_set = true;
}

void LQRAlgorithm::getSolution(Model::input_vector_t &u)
{
    assert(this->u);
    u = this->u.value();
}

void LQRAlgorithm::loadParameters()
{
    ParameterServer param(model->getParameterFolder() + "/LQR.info");

    Model::state_vector_t q;
    Model::input_vector_t r;
    param.loadMatrix("state_weights", q);
    param.loadMatrix("input_weights", r);
    setStateWeights(q);
    setInputWeights(r);
}

} // namespace scpp