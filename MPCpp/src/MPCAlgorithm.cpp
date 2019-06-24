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
    U.resize(Model::input_dim, K);

    Model::state_vector_t x_eq;
    Model::input_vector_t u_eq;
    model->getOperatingPoint(x_eq, u_eq);
    Eigen::MatrixXd X_eq(size_t(Model::state_dim), K);
    Eigen::MatrixXd U_eq(size_t(Model::input_dim), K);
    X_eq.col(0) = x_eq;
    X_eq.col(1) = x_eq;
    U_eq.col(0) = u_eq;
    U_eq.col(1) = u_eq;
    Model::state_matrix_v_t A_eq(K);
    Model::control_matrix_v_t B_eq(K);
    Model::control_matrix_v_t C_eq(K);
    Model::state_vector_v_t z_eq(K);

    double T;
    model->getTimeHorizon(T);

    multipleShooting(*model, T, X_eq, U_eq, A_eq, B_eq, C_eq, z_eq);
    A = A_eq.at(0);
    B = B_eq.at(0);
    C = C_eq.at(0);
    z = z_eq.at(0);

    socp = mpc::buildSCOP(*model, X, U, u_init, x_init, x_des, A, B, C, z);

    cacheIndices();

    solver = std::make_unique<EcosWrapper>(socp);
}

void MPCAlgorithm::setInitialInput(const Model::input_vector_t &u)
{
    u_init << u;
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
    // print("Solving model {}\n", Model::getModelName());

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
    for (size_t k = 0; k < K; k++)
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
    for (size_t k = 0; k < K; k++)
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
