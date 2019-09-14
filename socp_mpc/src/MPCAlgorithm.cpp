#include "MPCAlgorithm.hpp"
#include "timing.hpp"

using fmt::format;
using fmt::print;
using std::string;
using std::vector;

namespace mpc
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
    param.loadScalar("time_horizon", time_horizon);
    param.loadMatrix("state_weights_intermediate", state_weights_intermediate);
    param.loadMatrix("state_weights_terminal", state_weights_terminal);
    param.loadMatrix("input_weights", input_weights);

    setStateWeights(state_weights_intermediate, state_weights_terminal);
    setInputWeights(input_weights);
}

void MPCAlgorithm::initialize(bool constant_dynamics,
                              bool intermediate_cost_active)
{
    print("Initializing model '{}'.\n", Model::getModelName());

    X.resize(K);
    U.resize(K - 1);

    print("Computing dynamics.\n");
    const double timer_dynamics = tic();
    model->initializeModel();
    print("{:<{}}{:.2f}ms\n", "Time, dynamics:", 50, toc(timer_dynamics));

    print("Discretizing.\n");
    const double timer_discretize = tic();
    Model::state_vector_t x_eq;
    Model::input_vector_t u_eq;
    model->getOperatingPoint(x_eq, u_eq);
    const double dt = time_horizon / (K - 1);
    discretization::exactLinearDiscretization(model, dt, x_eq, u_eq, A, B, z);
    print("{:<{}}{:.2f}ms\n", "Time, discretization:", 50, toc(timer_discretize));

    // to check
    // Model::dynamic_matrix_t X_eq(10, K);
    // X_eq.col(0) = x_eq;
    // X_eq.col(1) = x_eq;
    // Model::dynamic_matrix_t U_eq(3, K);
    // U_eq.col(0) = u_eq;
    // U_eq.col(1) = u_eq;
    // Model::state_matrix_v_t A_bar(K);
    // Model::control_matrix_v_t B_bar(K);
    // Model::control_matrix_v_t C_bar(K);
    // Model::state_vector_v_t z_bar(K);

    // discretization::multipleShooting(model, time_horizon, X_eq, U_eq, A_bar, B_bar, C_bar, z_bar);
    // std::cout << (A_bar[0] - A).cwiseAbs().maxCoeff() << "\n";
    // std::cout << (B_bar[0] + C_bar[0] - B).cwiseAbs().maxCoeff() << "\n";
    // std::cout << (z_bar[0] - z).cwiseAbs().maxCoeff() << "\n";

    socp = mpc::buildSCOP(model,
                          X, U,
                          x_init, x_final,
                          state_weights_intermediate, state_weights_terminal, input_weights,
                          A, B, z,
                          constant_dynamics, intermediate_cost_active);

    cacheIndices();

    solver = std::make_unique<EcosWrapper>(socp);
}

void MPCAlgorithm::setInitialState(const Model::state_vector_t &x) { x_init = x; }

void MPCAlgorithm::setFinalState(const Model::state_vector_t &x) { x_final = x; }

void MPCAlgorithm::setStateWeights(const Model::state_vector_t &intermediate,
                                   const Model::state_vector_t &terminal)
{
    state_weights_intermediate = intermediate;
    state_weights_terminal = terminal;
}

void MPCAlgorithm::setInputWeights(const Model::input_vector_t &intermediate)
{
    input_weights = intermediate;
}

void MPCAlgorithm::solve()
{
    print("Solving model {}\n", Model::getModelName());
    const double timer_solve = tic();
    if (nondimensionalize)
    {
        model->nondimensionalize();
    }
    const int exit_flag = solver->solveProblem(false);
    if (exit_flag == -4)
    {
        print("Process interrupted.");
        std::terminate();
    }

    if (exit_flag != 0)
    {
        print("There was an issue (Exit code {}) while solving the problem. Continuing.\n", exit_flag);
    }

    if (nondimensionalize)
    {
        model->redimensionalize();
    }
    print("{:<{}}{:.2f}ms\n", "Time, solve:", 50, toc(timer_solve));

    readSolution();
}

void MPCAlgorithm::cacheIndices()
{
    // cache indices for performance
    X_indices.resize(Model::state_dim, K);
    U_indices.resize(Model::input_dim, K);
    for (size_t k = 0; k < K; k++)
    {
        for (size_t i = 0; i < Model::state_dim; i++)
        {
            X_indices(i, k) = socp.getTensorVariableIndex("X", {i, k});
        }
    }
    for (size_t k = 0; k < K - 1; k++)
    {
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            U_indices(i, k) = socp.getTensorVariableIndex("U", {i, k});
        }
    }
}

void MPCAlgorithm::readSolution()
{
    for (size_t k = 0; k < K; k++)
    {
        for (size_t i = 0; i < Model::state_dim; i++)
        {
            X.at(k)(i) = solver->getSolutionValue(X_indices(i, k));
        }
    }
    for (size_t k = 0; k < K - 1; k++)
    {
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            U.at(k)(i) = solver->getSolutionValue(U_indices(i, k));
        }
    }
}

void MPCAlgorithm::getSolution(Model::state_vector_v_t &X, Model::input_vector_v_t &U)
{
    X = this->X;
    U = this->U;
}

} // namespace mpc