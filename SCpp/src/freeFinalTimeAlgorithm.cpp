#include "freeFinalTimeAlgorithm.hpp"
#include "timing.hpp"

using fmt::format;
using fmt::print;
using std::string;
using std::vector;

FreeFinalTimeAlgorithm::FreeFinalTimeAlgorithm(std::shared_ptr<Model> model)
    : param("../SCpp/config/SCParameters.info"), model(model)
{
    loadParameters();
}

void FreeFinalTimeAlgorithm::loadParameters()
{
    param.loadScalar("K", K);

    param.loadScalar("delta_tol", delta_tol);
    param.loadScalar("max_iterations", max_iterations);
    param.loadScalar("nu_tol", nu_tol);

    param.loadScalar("trust_region_factor", trust_region_factor);
    param.loadScalar("weight_trust_region_time", weight_trust_region_time);
    param.loadScalar("weight_trust_region_trajectory", weight_trust_region_trajectory);
    param.loadScalar("weight_virtual_control", weight_virtual_control);
}

void FreeFinalTimeAlgorithm::initialize()
{
    model->initializeModel();

    A_bar.resize(K - 1);
    B_bar.resize(K - 1);
    C_bar.resize(K - 1);
    S_bar.resize(K - 1);
    z_bar.resize(K - 1);

    X.resize(Model::state_dim, K);
    U.resize(Model::input_dim, K);

    socp = sc::buildSCOP(*model,
                         weight_trust_region_time, weight_trust_region_trajectory, weight_virtual_control,
                         X, U, sigma, A_bar, B_bar, C_bar, S_bar, z_bar);

    cacheIndices();

    solver = std::make_unique<EcosWrapper>(socp);
}

bool FreeFinalTimeAlgorithm::iterate()
{
    // discretize
    const double timer_iteration = tic();
    double timer = tic();
    calculateDiscretization(*model, sigma, X, U, A_bar, B_bar, C_bar, S_bar, z_bar);
    print("{:<{}}{:.2f}ms\n", "Time, discretization:", 50, toc(timer));

    // solve the problem
    timer = tic();
    solver->solveProblem();
    print("{:<{}}{:.2f}ms\n", "Time, solver:", 50, toc(timer));

    readSolution();

    // // check feasibility
    timer = tic();
    if (!socp.feasibilityCheck(solver->getSolutionVector()))
    {
        throw std::runtime_error("Solver produced an invalid solution.\n");
    }
    print("{:<{}}{:.2f}ms\n", "Time, solution check:", 50, toc(timer));

    // print iteration summary
    print("\n");
    print("{:<{}}{: .4f}\n", "Trajectory Time", 50, sigma);
    print("{:<{}}{: .4f}\n", "Norm Virtual Control", 50, solver->getSolutionValue("norm1_nu", {}));
    print("{:<{}}{: .4f}\n", "State Input Delta", 50, solver->getSolutionValue("Delta_sigma", {}));
    print("{:<{}}{: .4f}\n\n", "Trust Region Delta", 50, solver->getSolutionValue("norm2_Delta", {}));

    print("{:<{}}{:.2f}ms\n\n", "Time, iteration:", 50, toc(timer_iteration));

    // check for convergence
    return solver->getSolutionValue("norm2_Delta", {}) < delta_tol && solver->getSolutionValue("norm1_nu", {}) < nu_tol;
}

void FreeFinalTimeAlgorithm::solve(bool warm_start)
{
    print("Solving model {}\n", Model::getModelName());

    model->nondimensionalize();

    if (warm_start)
    {
        model->nondimensionalizeTrajectory(X, U);
    }
    else
    {
        loadParameters();
        model->getInitializedTrajectory(X, U);
        model->getFinalTimeGuess(sigma);
    }

    Model::param_vector_t model_params;
    model->getNewModelParameters(model_params);
    model->updateParameters(model_params);

    const double timer_total = tic();

    size_t iteration = 0;
    bool converged = false;
    while (iteration <= max_iterations)
    {
        iteration++;
        print("{:=^{}}\n", format("<Iteration {}>", iteration), 60);

        converged = iterate();

        if (converged)
        {
            print("Converged after {} iterations.\n\n", iteration);
            break;
        }
        else if (iteration > 5)
        {
            // else increase trust region weight
            weight_trust_region_time *= trust_region_factor;
            weight_trust_region_trajectory *= trust_region_factor;
        }
    }

    if (not converged)
    {
        print("No convergence after {} iterations.\n\n", max_iterations);
    }

    model->redimensionalize();
    model->redimensionalizeTrajectory(X, U);

    print("{:<{}}{:.2f}ms\n", "Time, total:", 50, toc(timer_total));
}

void FreeFinalTimeAlgorithm::cacheIndices()
{
    // cache indices for performance
    sigma_index = socp.getTensorVariableIndex("sigma", {});
    X_indices.resizeLike(X);
    U_indices.resizeLike(U);
    for (size_t k = 0; k < K; k++)
    {
        for (size_t i = 0; i < Model::state_dim; i++)
        {
            X_indices(i, k) = socp.getTensorVariableIndex("X", {i, k});
        }
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            U_indices(i, k) = socp.getTensorVariableIndex("U", {i, k});
        }
    }
}

void FreeFinalTimeAlgorithm::readSolution()
{
    for (size_t k = 0; k < K; k++)
    {
        for (size_t i = 0; i < Model::state_dim; i++)
        {
            X(i, k) = solver->getSolutionValue(X_indices(i, k));
        }
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            U(i, k) = solver->getSolutionValue(U_indices(i, k));
        }
    }
    sigma = solver->getSolutionValue(sigma_index);
}

void FreeFinalTimeAlgorithm::getSolution(Model::dynamic_matrix_t &X, Model::dynamic_matrix_t &U, double &t)
{
    X = this->X;
    U = this->U;
    t = this->sigma;
}
