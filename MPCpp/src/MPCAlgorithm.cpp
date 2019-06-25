#include "MPCAlgorithm.hpp"
#include "timing.hpp"

using fmt::format;
using fmt::print;
using std::string;
using std::vector;

MPCAlgorithm::MPCAlgorithm(std::shared_ptr<Model> model)
    : param("../MPCpp/config/MPCParameters.info"), model(model)
{
    loadParameters();
}

void MPCAlgorithm::loadParameters()
{
    param.loadScalar("K", K);
}

void MPCAlgorithm::initialize()
{
    model->initializeModel();

    X.resize(Model::state_dim, K);
    U.resize(Model::input_dim, K - 1);

    Model::state_vector_t x_eq;
    Model::input_vector_t u_eq;
    model->getOperatingPoint(x_eq, u_eq);

    double T;
    model->getTimeHorizon(T);

    const double dt = T / (K - 1);

    print("Discretizing model {}\n", Model::getModelName());
    const double timer_discretize = tic();
    exactLinearDiscretization(*model, dt, x_eq, u_eq, A, B, z);
    print("{:<{}}{:.2f}ms\n", "Time, discretization:", 50, toc(timer_discretize));

    socp = mpc::buildSCOP(*model, X, U, x_init, x_des, A, B, z);

    cacheIndices();

    solver = std::make_unique<EcosWrapper>(socp);
}

void MPCAlgorithm::getTimeSteps(size_t &K)
{
    K = this->K;
}

void MPCAlgorithm::setInitialState(const Model::state_vector_t &x)
{
    x_init << x;
}

void MPCAlgorithm::setDesiredState(const Model::state_vector_t &x)
{
    x_des << x;
}

void MPCAlgorithm::solve()
{
    print("Solving model {}\n", Model::getModelName());
    const double timer_solve = tic();
    solver->solveProblem(false);
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
