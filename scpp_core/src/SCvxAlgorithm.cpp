#include "SCvxAlgorithm.hpp"
#include "simulation.hpp"
#include "timing.hpp"

using fmt::format;
using fmt::print;
using std::string;
using std::vector;

namespace scpp
{

SCvxAlgorithm::SCvxAlgorithm(Model::ptr_t model)
{
    this->model = model;
    loadParameters();

    all_X.reserve(max_iterations);
    all_U.reserve(max_iterations);
}

void SCvxAlgorithm::loadParameters()
{
    ParameterServer param(model->getParameterFolder() + "/SCvx.info");

    param.loadScalar("K", K);

    param.loadScalar("nondimensionalize", nondimensionalize);

    param.loadScalar("max_iterations", max_iterations);
    param.loadScalar("alpha", alpha);
    param.loadScalar("beta", beta);
    param.loadScalar("rho_0", rho_0);
    param.loadScalar("rho_1", rho_1);
    param.loadScalar("rho_2", rho_2);

    param.loadScalar("change_threshold", change_threshold);
    param.loadScalar("weight_virtual_control", weight_virtual_control);
    param.loadScalar("trust_region", trust_region);

    param.loadScalar("interpolate_input", interpolate_input);
}

void SCvxAlgorithm::initialize()
{
    print("Computing dynamics.\n");
    const double timer_dynamics = tic();
    model->initializeModel();
    print("{:<{}}{:.2f}ms\n", "Time, dynamics:", 50, toc(timer_dynamics));

    A_bar.resize(K - 1);
    B_bar.resize(K - 1);
    if (interpolate_input)
    {
        C_bar.resize(K - 1);
    }
    z_bar.resize(K - 1);

    X.resize(K);
    U.resize(interpolate_input ? K : K - 1);

    socp = buildSCvxProblem(model,
                            trust_region, weight_virtual_control,
                            X, U, A_bar, B_bar, C_bar, z_bar);

    cacheIndices();

    solver = std::make_unique<EcosWrapper>(socp);
}

bool SCvxAlgorithm::iterate()
{
    // discretize
    const double timer_iteration = tic();
    double timer = tic();

    if (interpolate_input)
    {
        discretization::multipleShooting(model, sigma, X, U, A_bar, B_bar, C_bar, z_bar);
    }
    else
    {
        discretization::multipleShooting(model, sigma, X, U, A_bar, B_bar, z_bar);
    }

    print("{:<{}}{:.2f}ms\n", "Time, discretization:", 50, toc(timer));

    bool converged = false;

    while (true)
    {
        // solve problem
        timer = tic();

        const Model::state_vector_v_t old_X = X;
        const Model::input_vector_v_t old_U = U;

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

        // compare linear and nonlinear costs
        const double nonlinear_cost_dynamics = getNonlinearCost();
        const double linear_cost_dynamics = solver->getSolutionValue("norm1_nu", {});

        // TODO: Consider linearized model constraints
        const double nonlinear_cost_constraints = 0.;
        const double linear_cost_constraints = 0.;

        const double nonlinear_cost = nonlinear_cost_dynamics + nonlinear_cost_constraints; // J
        const double linear_cost = linear_cost_dynamics + linear_cost_constraints;          // L

        if (not last_nonlinear_cost)
        {
            last_nonlinear_cost = nonlinear_cost;
            break;
        }

        const double actual_change = last_nonlinear_cost.value() - nonlinear_cost; // delta_J
        const double predicted_change = last_nonlinear_cost.value() - linear_cost; // delta_L

        last_nonlinear_cost = nonlinear_cost;

        print("{:<{}}{:.5f}\n", "Actual change:", 50, actual_change);
        print("{:<{}}{:.5f}\n", "Predicted change:", 50, predicted_change);

        if (std::abs(predicted_change) < change_threshold)
        {
            converged = true;
            break;
        }

        const double rho = actual_change / predicted_change;
        if (rho < rho_0)
        {
            trust_region /= alpha;
            print("Trust region too large. Solving again with radius={}\n", trust_region);
            print("--------------------------------------------------\n");
            X = old_X;
            U = old_U;
        }
        else
        {
            print("Solution accepted.\n");

            if (rho < rho_1)
            {
                print("Decreasing radius.\n");
                trust_region /= alpha;
            }
            else if (rho >= rho_2)
            {
                print("Increasing radius.\n");
                trust_region *= beta;
            }
            break;
        }
    }

    // print iteration summary
    print("\n");
    print("{:<{}}{: .4f}s\n\n", "Trajectory Time", 50, sigma);

    print("{:<{}}{:.2f}ms\n\n", "Time, iteration:", 50, toc(timer_iteration));

    return converged;
}

void SCvxAlgorithm::solve(bool warm_start)
{
    const double timer_total = tic();

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
        model->getInitializedTrajectory(X, U, sigma);
    }

    model->updateModelParameters();

    size_t iteration = 0;

    all_X.push_back(X);
    all_U.push_back(U);

    bool converged = false;

    while (iteration < max_iterations and not converged)
    {
        iteration++;
        print("{:=^{}}\n", format("<Iteration {}>", iteration), 60);

        converged = iterate();

        all_X.push_back(X);
        all_U.push_back(U);
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

void SCvxAlgorithm::cacheIndices()
{
    // cache indices for performance
    X_indices.resize(Model::state_dim, X.size());
    U_indices.resize(Model::input_dim, U.size());
    for (size_t k = 0; k < X.size(); k++)
    {
        for (size_t i = 0; i < Model::state_dim; i++)
        {
            X_indices(i, k) = socp.getTensorVariableIndex("X", {i, k});
        }
    }
    for (size_t k = 0; k < U.size(); k++)
    {
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            U_indices(i, k) = socp.getTensorVariableIndex("U", {i, k});
        }
    }
}

void SCvxAlgorithm::readSolution()
{
    for (size_t k = 0; k < X.size(); k++)
    {
        for (size_t i = 0; i < Model::state_dim; i++)
        {
            X[k](i) = solver->getSolutionValue(X_indices(i, k));
        }
    }
    for (size_t k = 0; k < U.size(); k++)
    {
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            U[k](i) = solver->getSolutionValue(U_indices(i, k));
        }
    }
}

void SCvxAlgorithm::getSolution(Model::state_vector_v_t &X, Model::input_vector_v_t &U, double &t)
{
    X = this->X;
    U = this->U;
    t = sigma;
}

void SCvxAlgorithm::getAllSolutions(std::vector<Model::state_vector_v_t> &X,
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

    t.resize(X.size());
    std::fill(t.begin(), t.end(), sigma);
}

double SCvxAlgorithm::getNonlinearCost()
{
    double nonlinear_cost_dynamics = 0.;
    for (size_t k = 0; k < K - 1; k++)
    {
        Model::state_vector_t x = X.at(k);
        const Model::input_vector_t u0 = U.at(k);
        const Model::input_vector_t u1 = interpolate_input ? U.at(k + 1) : u0;

        simulate(model, sigma / (K - 1), u0, u1, x);

        const double virtual_control_cost = (x - X.at(k + 1)).lpNorm<1>();
        nonlinear_cost_dynamics += virtual_control_cost;
    }

    return nonlinear_cost_dynamics;
}

} // namespace scpp