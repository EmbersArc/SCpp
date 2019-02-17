#include "freeFinalTimeAlgorithm.hpp"

using fmt::format;
using fmt::print;
using std::ofstream;
using std::string;
using std::vector;

namespace fs = std::filesystem;

freeFinalTimeAlgorithm::freeFinalTimeAlgorithm(std::shared_ptr<Model> model)
    : param("../SCpp/config/SCParameters.info"), model(model)
{
}

void freeFinalTimeAlgorithm::loadParameters()
{
    param.loadScalar("K", K);

    param.loadScalar("weight_trust_region_time", weight_trust_region_time);
    param.loadScalar("weight_trust_region_trajectory", weight_trust_region_trajectory);
    param.loadScalar("weight_virtual_control", weight_virtual_control);
    param.loadScalar("trust_region_factor", trust_region_factor);
    param.loadScalar("nu_tol", nu_tol);
    param.loadScalar("delta_tol", delta_tol);
    param.loadScalar("max_iterations", max_iterations);
}

void freeFinalTimeAlgorithm::initialize()
{
    loadParameters();

    model->nondimensionalize();
    model->initializeModel();

    A_bar.resize(K - 1);
    B_bar.resize(K - 1);
    C_bar.resize(K - 1);
    S_bar.resize(K - 1);
    z_bar.resize(K - 1);

    X.resize(Model::state_dim, K);
    U.resize(Model::input_dim, K);
    model->initializeTrajectory(X, U);
    sigma = model->getFinalTimeGuess();

    Model::param_vector_t model_params;
    model->getNewModelParameters(model_params);
    model->updateParameters(model_params);

    socp = sc::buildSCOP(*model,
                         weight_trust_region_time, weight_trust_region_trajectory, weight_virtual_control,
                         X, U, sigma, A_bar, B_bar, C_bar, S_bar, z_bar);

    cacheIndices();

    solver = std::make_unique<EcosWrapper>(socp);
}

bool freeFinalTimeAlgorithm::iterate()
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

    getSolution();

    // print iteration summary
    print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
    print("\n");
    print("{:<{}}{: .4f}\n", "sigma", 50, sigma);
    print("{:<{}}{: .4f}\n", "norm1_nu", 50, solver->getSolutionValue("norm1_nu", {}));
    print("{:<{}}{: .4f}\n", "Delta_sigma", 50, solver->getSolutionValue("Delta_sigma", {}));
    print("{:<{}}{: .4f}\n", "norm2_Delta", 50, solver->getSolutionValue("norm2_Delta", {}));
    print("\n");
    print("{:<{}}{:.2f}ms\n", "Time, iteration:", 50, toc(timer_iteration));
    print("\n");

    // check for convergence
    return solver->getSolutionValue("norm2_Delta", {}) < delta_tol && solver->getSolutionValue("norm1_nu", {}) < nu_tol;
}

void freeFinalTimeAlgorithm::solve()
{
    const double timer_total = tic();

    size_t iteration = 0;
    bool converged = false;
    while (iteration <= max_iterations)
    {
        iteration++;
        string itString = format("<Iteration {}>", iteration);
        print("{:=^{}}\n", itString, 60);

        converged = iterate();
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
    else
    {
        // write solution to files
        double timer = tic();
        string outputDirectory = format("{}/{}", getOutputPath(), time(0));

        if (not fs::exists(outputDirectory) and not fs::create_directories(outputDirectory))
        {
            throw std::runtime_error("Could not create output directory!");
        }

        model->redimensionalizeTrajectory(X, U);
        {
            ofstream f(outputDirectory + "/X.txt");
            f << X;
        }
        {
            ofstream f(outputDirectory + "/U.txt");
            f << U;
        }
        {
            ofstream f(outputDirectory + "/t.txt");
            f << sigma;
        }
        print("{:<{}}{:.2f}ms\n", "Time, solution file:", 50, toc(timer));
    }
    print("{:<{}}{:.2f}ms\n", "Time, total:", 50, toc(timer_total));
}

string freeFinalTimeAlgorithm::getOutputPath()
{
    return format("../output/{}", Model::getModelName());
}

void freeFinalTimeAlgorithm::cacheIndices()
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

void freeFinalTimeAlgorithm::getSolution()
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