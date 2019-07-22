#include "MPCAlgorithm.hpp"
#include "timing.hpp"

using fmt::format;
using fmt::print;
using std::string;
using std::vector;

namespace mpc
{

MPCAlgorithm::MPCAlgorithm(std::shared_ptr<Model> model, const std::string parameter_path) : model(model)
{
    std::string filename = "MPCParameters.info";
    std::string path;
    if (parameter_path == "")
    {
        path = "../socp_mpc/config/" + filename;
    }
    else
    {
        path = parameter_path + filename;
    }
    loadParameters(path);
}

void MPCAlgorithm::loadParameters(const std::string &path)
{
    ParameterServer param(path);

    param.loadScalar("K", K);
    param.loadScalar("nondimensionalize", nondimensionalize);
}

void MPCAlgorithm::initialize(bool constant_dynamics)
{
    print("Initializing model '{}'.\n", Model::getModelName());

    X.resize(Model::state_dim, K);
    U.resize(Model::input_dim, K - 1);

    model->getStateWeights(state_weights_intermediate, state_weights_terminal);
    model->getInputWeights(input_weights);

    double T;
    model->getTimeHorizon(T);

    dt = T / (K - 1);

    print("Computing dynamics.\n");
    const double timer_dynamics = tic();
    model->initializeModel();
    print("{:<{}}{:.2f}ms\n", "Time, dynamics:", 50, toc(timer_dynamics));

    Model::state_vector_t x_eq;
    Model::input_vector_t u_eq;
    model->getOperatingPoint(x_eq, u_eq);
    discretize(x_eq, u_eq);

    socp = mpc::buildSCOP(*model,
                          X, U,
                          x_init, x_final,
                          state_weights_intermediate, state_weights_terminal, input_weights,
                          A, B, z,
                          constant_dynamics);

    cacheIndices();

    solver = std::make_unique<EcosWrapper>(socp);
}

void MPCAlgorithm::discretize(const Model::state_vector_t &x, const Model::input_vector_t &u)
{
    print("Discretizing.\n");
    const double timer_discretize = tic();
    exactLinearDiscretization(*model, dt, x, u, A, B, z);
    print("{:<{}}{:.2f}ms\n", "Time, discretization:", 50, toc(timer_discretize));
}

void MPCAlgorithm::getTimeSteps(size_t &K) { K = this->K; }

void MPCAlgorithm::setInitialState(const Model::state_vector_t &x) { x_init << x; }

void MPCAlgorithm::setFinalState(const Model::state_vector_t &x) { x_final << x; }

void MPCAlgorithm::setStateWeights(const Model::state_vector_t &intermediate, const Model::state_vector_t &terminal)
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
    solver->solveProblem(false);
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
    X_indices.resizeLike(X);
    U_indices.resizeLike(U);
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
            X(i, k) = solver->getSolutionValue(X_indices(i, k));
        }
    }
    for (size_t k = 0; k < K - 1; k++)
    {
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            U(i, k) = solver->getSolutionValue(U_indices(i, k));
        }
    }
}

void MPCAlgorithm::getSolution(Model::dynamic_matrix_t &X, Model::dynamic_matrix_t &U)
{
    X = this->X;
    U = this->U;
}

} // namespace mpc