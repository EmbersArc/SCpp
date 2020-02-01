#include "SCAlgorithm.hpp"
#include "timing.hpp"
#include "simulation.hpp"
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
    solver = std::make_unique<op::Solver>(socp);
}

bool SCAlgorithm::iterate()
{
    // discretize
    const double timer_iteration = tic();
    double timer = tic();
    discretization::multipleShooting(model, td, dd);
    print("{:<{}}{:.2f}ms\n", "Time, discretization:", 50, toc(timer));

    // solve the problem
    print("\n");
    print("Solving problem.\n");
    timer = tic();
    const bool success = solver->solveProblem(false);
    print("Solver message:\n");
    print("> {}\n", solver->getResultString());
    print("{:<{}}{:.2f}ms\n", "Time, solver:", 50, toc(timer));
    print("\n");

    timer = tic();
    std::vector<bool> defects = calculateDefects();
    print("{:<{}}{:.2f}ms\n", "Time, defects:", 50, toc(timer));
    print("Defect pattern:\n");
    for (bool defect : defects)
    {
        print("{}", defect ? "x" : "-");
    }
    print("\n\n");

    if (not success)
    {
        print("Solver failed to find a solution. Terminating.\n");
        std::terminate();
    }

    readSolution();

    // check feasibility
    timer = tic();
    if (!socp.isFeasible())
    {
        throw std::runtime_error("Solver produced an invalid solution.\n");
    }
    print("{:<{}}{:.2f}ms\n", "Time, solution check:", 50, toc(timer));

    double norm1_nu, norm1_delta;
    socp.readSolution("norm1_nu", norm1_nu);
    socp.readSolution("norm1_delta", norm1_delta);

    if (norm1_nu < nu_tol)
    {
        weight_trust_region_trajectory *= 2.;
    }

    // print iteration summary
    print("\n");
    print("{:<{}}{: .4f}\n", "Norm Virtual Control", 50, norm1_nu);
    if (free_final_time)
    {
        double delta_sigma;
        socp.readSolution("delta_sigma", delta_sigma);
        print("{:<{}}{: .4f}\n", "Time Trust Region Delta", 50, delta_sigma);
    }
    print("{:<{}}{: .4f}\n\n", "Trajectory Trust Region Delta", 50, norm1_delta);

    print("{:<{}}{: .4f}s\n\n", "Trajectory Time", 50, td.t);

    print("{:<{}}{:.2f}ms\n\n", "Time, iteration:", 50, toc(timer_iteration));

    // check for convergence
    return norm1_delta < delta_tol and norm1_nu < nu_tol;
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

    while (iteration < max_iterations and not converged)
    {
        iteration++;

        print("{:=^{}}\n", format("<Iteration {}>", iteration), 60);
        converged = iterate();

        all_td.push_back(td);
    }

    print("{:=^{}}\n\n", "", 60);

    if (converged)
    {
        print("Converged after {} iterations.\n\n", iteration);
    }
    else
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

void SCAlgorithm::readSolution()
{
    if (free_final_time)
    {
        socp.readSolution("sigma", td.t);
    }
    Eigen::MatrixXd X, U;
    socp.readSolution("X", X);
    socp.readSolution("U", U);

    for (size_t k = 0; k < td.n_X(); k++)
    {
        td.X[k] = X.col(k);
    }
    for (size_t k = 0; k < td.n_U(); k++)
    {
        td.U[k] = U.col(k);
    }
}

void SCAlgorithm::getSolution(trajectory_data_t &trajectory) const
{
    trajectory = td;
}

void SCAlgorithm::getAllSolutions(std::vector<trajectory_data_t> &all_trajectories)
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

std::vector<bool> SCAlgorithm::calculateDefects()
{
    std::vector<bool> pattern;

    for (size_t k = 0; k < K - 1; k++)
    {
        Model::state_vector_t x = td.X.at(k);
        const Model::input_vector_t u0 = td.U.at(k);
        const Model::input_vector_t u1 = interpolate_input ? td.U.at(k + 1) : u0;

        const double dt = td.t / (K - 1);
        simulate(model, dt, u0, u1, x);

        const double defect = (x - td.X.at(k + 1)).squaredNorm();

        pattern.push_back(defect > nu_tol);
    }

    return pattern;
}

} // namespace scpp