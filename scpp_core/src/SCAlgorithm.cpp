#include "SCAlgorithm.hpp"
#include "timing.hpp"
#include "discretization.hpp"

using fmt::format;
using fmt::print;
using std::string;
using std::vector;

namespace scpp
{

SCAlgorithm::SCAlgorithm(Model::ptr_t model)
{
    this->model = model;
    loadParameters();

    all_td.reserve(max_iterations);
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

    param.loadScalar("weight_time", weight_time);
    param.loadScalar("weight_virtual_control", weight_virtual_control);
    param.loadScalar("trust_region_factor", trust_region_factor);
    param.loadScalar("weight_trust_region_trajectory", weight_trust_region_trajectory);

    param.loadScalar("interpolate_input", interpolate_input);

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

    dd.initialize(K, interpolate_input, free_final_time);

    td.initialize(K, interpolate_input);

    socp = buildSCProblem(weight_time, weight_trust_region_time,
                          weight_trust_region_trajectory, weight_virtual_control,
                          td, dd);
    model->addApplicationConstraints(socp, td.X, td.U);
    cacheIndices();

    solver = std::make_unique<EcosWrapper>(socp);
}

bool SCAlgorithm::iterate()
{
    // discretize
    const double timer_iteration = tic();
    double timer = tic();
    discretization::multipleShooting(model, td, dd);

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

    print("{:<{}}{: .4f}s\n\n", "Trajectory Time", 50, td.t);

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
            model->nondimensionalizeTrajectory(td);
    }
    else
    {
        loadParameters();
        model->getInitializedTrajectory(td);
    }

    model->updateModelParameters();

    const double timer_total = tic();

    size_t iteration = 0;
    bool converged = false;

    all_td.push_back(td);

    while (iteration < max_iterations)
    {
        iteration++;
        print("{:=^{}}\n", format("<Iteration {}>", iteration), 60);

        converged = iterate();

        all_td.push_back(td);

        if (converged)
        {
            print("Converged after {} iterations.\n\n", iteration);
            break;
        }
        else if (iteration > 2)
        {
            // else increase trust region weight
            // weight_trust_region_time *= trust_region_factor;
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
        model->redimensionalizeTrajectory(td);
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
    X_indices.resize(Model::state_dim, td.n_X());
    U_indices.resize(Model::input_dim, td.n_U());
    for (size_t k = 0; k < td.n_X(); k++)
    {
        for (size_t i = 0; i < Model::state_dim; i++)
        {
            X_indices(i, k) = socp.getTensorVariableIndex("X", {i, k});
        }
    }
    for (size_t k = 0; k < td.n_U(); k++)
    {
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
        td.t = solver->getSolutionValue(sigma_index);
    }
    for (size_t k = 0; k < td.n_X(); k++)
    {
        for (size_t i = 0; i < Model::state_dim; i++)
        {
            td.X[k](i) = solver->getSolutionValue(X_indices(i, k));
        }
    }
    for (size_t k = 0; k < td.n_U(); k++)
    {
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            td.U[k](i) = solver->getSolutionValue(U_indices(i, k));
        }
    }
}

void SCAlgorithm::getSolution(TrajectoryData &trajectory) const
{
    trajectory = td;
}

void SCAlgorithm::getAllSolutions(std::vector<TrajectoryData> &all_trajectories)
{
    if (nondimensionalize)
    {
        auto all_td_redim = all_td;
        for (auto &td_redim : all_td_redim)
        {
            model->redimensionalizeTrajectory(td_redim);
        }
        all_trajectories = all_td_redim;
    }
    else
    {
        all_trajectories = all_td;
    }
}

} // namespace scpp