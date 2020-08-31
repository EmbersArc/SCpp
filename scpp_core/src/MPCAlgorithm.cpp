#include "MPCAlgorithm.hpp"
#include "MPCProblem.hpp"
#include "timing.hpp"
#include "discretization.hpp"

using fmt::print;

namespace scpp
{

MPCAlgorithm::MPCAlgorithm(Model::ptr_t model)
{
    this->model = model;
    loadParameters();
}

void MPCAlgorithm::loadParameters()
{
    ParameterServer param(model->getParameterFolder() + "/MPC.info");

    param.loadScalar("K", K);
    param.loadScalar("nondimensionalize", nondimensionalize);
    param.loadScalar("constant_dynamics", constant_dynamics);
    param.loadScalar("intermediate_cost_active", intermediate_cost_active);
    param.loadScalar("time_horizon", time_horizon);
    param.loadMatrix("state_weights_intermediate", state_weights_intermediate);
    param.loadMatrix("state_weights_terminal", state_weights_terminal);
    param.loadMatrix("input_weights", input_weights);

    setStateWeights(state_weights_intermediate, state_weights_terminal);
    setInputWeights(input_weights);
}

void MPCAlgorithm::initialize()
{
    print("[MPC] Starting controller for model '{}'.\n", Model::getModelName());

    assert(state_weights_set and input_weights_set);

    X.resize(K);
    U.resize(K - 1);

    model->updateModelParameters();

    print("[MPC] Discretizing.\n");
    const double timer_discretize = tic();
    Model::state_vector_t x_eq;
    Model::input_vector_t u_eq;
    model->getOperatingPoint(x_eq, u_eq);
    const double dt = time_horizon / (K - 1);
    scpp::discretization::exactLinearDiscretization(model, dt, x_eq, u_eq, A, B, z);
    print("{:<{}}{:.2f}ms\n", "[MPC] Time, discretization:", 50, toc(timer_discretize));

    socp = buildMPCProblem(X, U,
                           x_init, x_final,
                           state_weights_intermediate, state_weights_terminal, input_weights,
                           A, B, z,
                           constant_dynamics, intermediate_cost_active);
    model->addApplicationConstraints(socp, X, U);

    solver = std::make_unique<cvx::ecos::ECOSSolver>(*socp);

    initialized = true;
    print("[MPC] Controller started.\n");
}

void MPCAlgorithm::setInitialState(const Model::state_vector_t &x)
{
    x_init = x;
}

void MPCAlgorithm::setFinalState(const Model::state_vector_t &x)
{
    x_final = x;
}

void MPCAlgorithm::setStateWeights(const Model::state_vector_t &intermediate,
                                   const Model::state_vector_t &terminal)
{
    state_weights_intermediate = intermediate;
    state_weights_terminal = terminal;

    state_weights_set = true;
}

void MPCAlgorithm::setInputWeights(const Model::input_vector_t &intermediate)
{
    input_weights = intermediate;

    input_weights_set = true;
}

void MPCAlgorithm::solve()
{
    assert(initialized);

    print("Solving problem.\n");
    const double timer_solve = tic();
    if (nondimensionalize)
    {
        model->nondimensionalize();
    }

    solver->solve(false);

    print("Solver message:\n");
    print("> {}\n", solver->getResultString());

    if (nondimensionalize)
    {
        model->redimensionalize();
    }
    print("{:<{}}{:.2f}ms\n", "Time, solve:", 50, toc(timer_solve));

    readSolution();
}

void MPCAlgorithm::readSolution()
{
    cvx::MatrixX v_X, v_U;
    socp->getVariable("X", v_X);
    socp->getVariable("U", v_U);
    Eigen::MatrixXd X = eval(v_X);
    Eigen::MatrixXd U = eval(v_U);

    for (size_t k = 0; k < K; k++)
    {
        this->X[k] = X.col(k);
    }
    for (size_t k = 0; k < K-1; k++)
    {
        this->U[k] = U.col(k);
    }
}

void MPCAlgorithm::getSolution(Model::state_vector_v_t &X, Model::input_vector_v_t &U)
{
    X = this->X;
    U = this->U;
}

} // namespace scpp