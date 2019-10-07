#include "SCAlgorithm.hpp"
#include "timing.hpp"

using fmt::format;
using fmt::print;
using std::string;
using std::vector;

SCAlgorithm::SCAlgorithm(Model::ptr_t model)
{
    this->model = model;
    loadParameters();

    all_X.reserve(max_iterations);
    all_U.reserve(max_iterations);
    all_times.reserve(max_iterations);
}

void SCAlgorithm::loadParameters()
{
    ParameterServer param(model->getParameterFolder() + "/SC.info");

    param.loadScalar("K", K);

    param.loadScalar("free_final_time", free_final_time);

    param.loadScalar("nondimensionalize", nondimensionalize);

    param.loadScalar("delta_tol", delta_tol);
    param.loadScalar("max_iterations", max_iterations);
    param.loadScalar("nu_tol", nu_tol);

    param.loadScalar("weight_virtual_control", weight_virtual_control);
    param.loadScalar("trust_region_factor", trust_region_factor);
    param.loadScalar("weight_trust_region_trajectory", weight_trust_region_trajectory);

    if (free_final_time)
    {
        param.loadScalar("weight_trust_region_time", weight_trust_region_time);
    }
}

void SCAlgorithm::initialize()
{
    print("Computing dynamics.\n");
    const double timer_dynamics = tic();
    model->initializeModel();
    print("{:<{}}{:.2f}ms\n", "Time, dynamics:", 50, toc(timer_dynamics));

    A_bar.resize(K - 1);
    B_bar.resize(K - 1);
    C_bar.resize(K - 1);
    S_bar.resize(K - 1);
    z_bar.resize(K - 1);

    X.resize(K);
    U.resize(K);

    socp = sc::buildSCOP(model,
                         weight_trust_region_time, weight_trust_region_trajectory, weight_virtual_control,
                         X, U, sigma, A_bar, B_bar, C_bar, S_bar, z_bar, free_final_time);

    cacheIndices();

    solver = std::make_unique<EcosWrapper>(socp);
}

bool SCAlgorithm::iterate()
{
    // discretize
    const double timer_iteration = tic();
    double timer = tic();
    if (free_final_time)
    {
        discretization::multipleShootingVariableTime(model, sigma, X, U,
                                                     A_bar, B_bar, C_bar, S_bar, z_bar);
    }
    else
    {
        discretization::multipleShooting(model, sigma, X, U,
                                         A_bar, B_bar, C_bar, z_bar);
    }

    print("{:<{}}{:.2f}ms\n", "Time, discretization:", 50, toc(timer));

    // solve the problem
    timer = tic();
    solver->solveProblem(false);
    print("{:<{}}{:.2f}ms\n", "Time, solver:", 50, toc(timer));

    readSolution();

    // check feasibility
    timer = tic();
    if (!socp.feasibilityCheck(solver->getSolutionVector()))
    {
        throw std::runtime_error("Solver produced an invalid solution.\n");
    }
    print("{:<{}}{:.2f}ms\n", "Time, solution check:", 50, toc(timer));

    // print iteration summary
    print("\n");
    print("{:<{}}{: .4f}\n", "Norm Virtual Control", 50, solver->getSolutionValue("norm1_nu", {}));
    if (free_final_time)
    {
        print("{:<{}}{: .4f}\n", "State Input Delta", 50, solver->getSolutionValue("Delta_sigma", {}));
    }
    print("{:<{}}{: .4f}\n\n", "Trust Region Delta", 50, solver->getSolutionValue("norm2_Delta", {}));

    print("{:<{}}{: .4f}s\n\n", "Trajectory Time", 50, sigma);

    print("{:<{}}{:.2f}ms\n\n", "Time, iteration:", 50, toc(timer_iteration));

    // check for convergence
    return solver->getSolutionValue("norm2_Delta", {}) < delta_tol && solver->getSolutionValue("norm1_nu", {}) < nu_tol;
}

void SCAlgorithm::solve(bool warm_start)
{
    print("Solving model {}\n", Model::getModelName());

    if (nondimensionalize)
        model->nondimensionalize();

    if (warm_start)
    {
        if (nondimensionalize)
            model->nondimensionalizeTrajectory(X, U);
    }
    else
    {
        loadParameters();
        model->getInitializedTrajectory(X, U);
    }

    Model::param_vector_t model_params;
    model->getNewModelParameters(model_params);
    model->updateParameters(model_params);

    const double timer_total = tic();

    size_t iteration = 0;
    bool converged = false;

    all_X.push_back(X);
    all_U.push_back(U);
    all_times.push_back(sigma);

    while (iteration < max_iterations)
    {
        iteration++;
        print("{:=^{}}\n", format("<Iteration {}>", iteration), 60);

        converged = iterate();

        all_X.push_back(X);
        all_U.push_back(U);
        all_times.push_back(sigma);

        if (converged)
        {
            print("Converged after {} iterations.\n\n", iteration);
            break;
        }
        else if (iteration > 2)
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

    if (nondimensionalize)
    {
        model->redimensionalize();
        model->redimensionalizeTrajectory(X, U);
    }
    print("{:<{}}{:.2f}ms\n", "Time, total:", 50, toc(timer_total));
}

void SCAlgorithm::cacheIndices()
{
    // cache indices for performance
    if (free_final_time)
    {
        sigma_index = socp.getTensorVariableIndex("sigma", {});
    }
    X_indices.resize(Model::state_dim, X.size());
    U_indices.resize(Model::input_dim, X.size());
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

void SCAlgorithm::readSolution()
{
    if (free_final_time)
    {
        sigma = solver->getSolutionValue(sigma_index);
    }
    for (size_t k = 0; k < K; k++)
    {
        for (size_t i = 0; i < Model::state_dim; i++)
        {
            X[k](i) = solver->getSolutionValue(X_indices(i, k));
        }
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            U[k](i) = solver->getSolutionValue(U_indices(i, k));
        }
    }
}

void SCAlgorithm::getSolution(Model::state_vector_v_t &X, Model::input_vector_v_t &U, double &t)
{
    X = this->X;
    U = this->U;
    t = this->sigma;
}

void SCAlgorithm::getAllSolutions(std::vector<Model::state_vector_v_t> &X,
                                  std::vector<Model::input_vector_v_t> &U,
                                  std::vector<double> &t)
{
    if (nondimensionalize)
    {
        auto all_X_redim = all_X;
        auto all_U_redim = all_U;
        for (size_t k = 0; k < all_X.size(); k++)
        {
            model->redimensionalizeTrajectory(all_X_redim.at(k), all_U_redim.at(k));
        }
        X = all_X_redim;
        U = all_U_redim;
    }
    else
    {
        X = all_X;
        U = all_U;
    }

    t = all_times;
}
