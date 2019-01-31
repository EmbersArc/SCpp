/*
 * 
 * Independent implementation of 
 *     Successive Convexification for 6-DoF Mars Rocket Powered Landing 
 *     with Free-Final-Time (Michael Szmuk, Behcet Acikmese)
 * 
 * https://arxiv.org/abs/1802.03827
 * 
 */

#include <fstream>
#include <filesystem>

#include <fmt/format.h>

#include "activeModel.hpp"
#include "ecosWrapper.hpp"
#include "discretization.hpp"
#include "successiveConvexificationSOCP.hpp"
#include "timing.hpp"
#include "parameterServer.hpp"

using fmt::format;
using fmt::print;
using std::ofstream;
namespace fs = std::filesystem;

string getOutputPath()
{
    return format("../output/{}", Model::getModelName());
}

int main()
{
    ParameterServer param("../include/config/SCParameters.info");

    size_t K;
    param.loadScalar("K", K);

    print("Initializing model.\n");
    Model model;
    bool nondimensionalize;
    param.loadScalar("nondimensionalize", nondimensionalize);
    if (nondimensionalize)
    {
        model.nondimensionalize();
    }
    model.initializeModel();

    print("Initializing output directory.\n");
    if (fs::exists(getOutputPath()))
    {
        fs::remove_all(getOutputPath());
    }

    if (not fs::create_directories(getOutputPath()))
    {
        throw std::runtime_error("Could not create output directory!");
    }

    print("Initializing algorithm.\n");
    double weight_trust_region_time;
    double weight_trust_region_trajectory;
    double weight_virtual_control;
    double trust_region_factor;
    double nu_tol;
    double delta_tol;
    size_t max_iterations;

    param.loadScalar("weight_trust_region_time", weight_trust_region_time);
    param.loadScalar("weight_trust_region_trajectory", weight_trust_region_trajectory);
    param.loadScalar("weight_virtual_control", weight_virtual_control);
    param.loadScalar("trust_region_factor", trust_region_factor);
    param.loadScalar("nu_tol", nu_tol);
    param.loadScalar("delta_tol", delta_tol);
    param.loadScalar("max_iterations", max_iterations);

    Model::state_matrix_v_t A_bar(K - 1);
    Model::control_matrix_v_t B_bar(K - 1);
    Model::control_matrix_v_t C_bar(K - 1);
    Model::state_vector_v_t S_bar(K - 1);
    Model::state_vector_v_t z_bar(K - 1);

    print("Initializing trajectory.\n");
    Eigen::MatrixXd X(size_t(Model::state_dim_), K);
    Eigen::MatrixXd U(size_t(Model::input_dim_), K);
    model.initializeTrajectory(X, U);
    double sigma = model.getFinalTimeGuess();

    print("Initializing solver.\n");
    op::SecondOrderConeProgram socp = sc::buildSCSOCP(model,
                                                      weight_trust_region_time, weight_trust_region_trajectory, weight_virtual_control,
                                                      X, U, sigma, A_bar, B_bar, C_bar, S_bar, z_bar);

    // Cache indices for performance
    const size_t sigma_index = socp.getTensorVariableIndex("sigma", {});
    Eigen::MatrixXi X_indices(size_t(Model::state_dim_), K);
    Eigen::MatrixXi U_indices(size_t(Model::input_dim_), K);
    for (size_t k = 0; k < K; k++)
    {
        for (size_t i = 0; i < Model::state_dim_; i++)
        {
            X_indices(i, k) = socp.getTensorVariableIndex("X", {i, k});
        }
        for (size_t i = 0; i < Model::input_dim_; i++)
        {
            U_indices(i, k) = socp.getTensorVariableIndex("U", {i, k});
        }
    }

    EcosWrapper solver(socp);

    print("Starting Successive Convexification.\n");
    const double timer_total = tic();
    for (size_t it = 0; it < max_iterations; it++)
    {
        string itString = format("<Iteration {}>", it);
        print("{:=^{}}\n", itString, 60);

        const double timer_iteration = tic();
        double timer = tic();
        calculateDiscretization(model, sigma, X, U, A_bar, B_bar, C_bar, S_bar, z_bar);
        print("{:<{}}{:.2f}ms\n", "Time, discretization:", 50, toc(timer));

        // Write solution to files
        timer = tic();
        Eigen::MatrixXd X_out(X);
        Eigen::MatrixXd U_out(U);
        model.redimensionalizeTrajectory(X_out, U_out);
        string file_name_prefix = format("{}/{:03}_", getOutputPath(), it);
        {
            ofstream f(file_name_prefix + "X.txt");
            f << X_out;
        }
        {
            ofstream f(file_name_prefix + "U.txt");
            f << U_out;
        }
        {
            ofstream f(file_name_prefix + "t.txt");
            f << sigma;
        }
        print("{:<{}}{:.2f}ms\n", "Time, solution file:", 50, toc(timer));

        // Solve the problem
        timer = tic();
        solver.solveProblem();
        print("{:<{}}{:.2f}ms\n", "Time, solver:", 50, toc(timer));

        // Check feasibility
        timer = tic();
        if (!socp.feasibilityCheck(solver.getSolutionVector()))
        {
            print("ERROR: Solver produced an invalid solution.\n");
            return EXIT_FAILURE;
        }
        print("{:<{}}{:.2f}ms\n", "Time, solution check:", 50, toc(timer));

        // Read solution
        sigma = solver.getSolutionValue(sigma_index);
        for (size_t k = 0; k < K; k++)
        {
            for (size_t i = 0; i < Model::state_dim_; i++)
            {
                X(i, k) = solver.getSolutionValue(X_indices(i, k));
            }
            for (size_t i = 0; i < Model::input_dim_; i++)
            {
                U(i, k) = solver.getSolutionValue(U_indices(i, k));
            }
        }

        // Output iteration summary
        print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
        print("\n");
        print("{:<{}}{: .4f}\n", "sigma", 50, sigma);
        print("{:<{}}{: .4f}\n", "norm1_nu", 50, solver.getSolutionValue("norm1_nu", {}));
        print("{:<{}}{: .4f}\n", "Delta_sigma", 50, solver.getSolutionValue("Delta_sigma", {}));
        print("{:<{}}{: .4f}\n", "norm2_Delta", 50, solver.getSolutionValue("norm2_Delta", {}));
        print("\n");
        print("{:<{}}{:.2f}ms\n", "Time, iteration:", 50, toc(timer_iteration));
        print("\n");
        if (solver.getSolutionValue("norm2_Delta", {}) < delta_tol && solver.getSolutionValue("norm1_nu", {}) < nu_tol)
        {
            print("Converged after {} iterations.\n", it + 1);
            break;
        }
        else if (solver.getSolutionValue("norm1_nu", {}) < nu_tol)
        {
            weight_trust_region_time *= trust_region_factor;
            weight_trust_region_trajectory *= trust_region_factor;
        }
        else if (it == max_iterations)
        {
            print("No convergence after {} iterations.\n", max_iterations);
        }
    }

    print("{:<{}}{:.2f}ms\n", "Time, total:", 50, toc(timer_total));
}
