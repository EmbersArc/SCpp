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

#include "active_model.hpp"
#include "ecosWrapper.hpp"
#include "discretization.hpp"
#include "successiveConvexificationSOCP.hpp"
#include "timing.hpp"

using fmt::format;
using fmt::print;
using std::array;
using std::ofstream;
using std::filesystem::create_directory;
using std::filesystem::remove_all;

string get_output_path()
{
    return format("../output/{}/", Model::getModelName());
}

int main()
{
    print("Initializing model.\n");
    Model model;
    model.initializeModel();

    print("Initializing path.\n");
    remove_all(get_output_path());
    create_directory(get_output_path());

    print("Initializing algorithm.\n");
    double weight_trust_region_sigma = 1e-1;
    double weight_trust_region_xu = 1e-3;
    double weight_virtual_control = 1e5;
    double nu_tol = 1e-8;
    double delta_tol = 1e-3;

    Eigen::Matrix<double, Model::state_dim_, K> X;
    Eigen::Matrix<double, Model::input_dim_, K> U;

    model.initializeTrajectory(X, U);

    double sigma = model.getFinalTimeGuess();

    array<Model::state_matrix_t, (K - 1)> A_bar;
    array<Model::control_matrix_t, (K - 1)> B_bar;
    array<Model::control_matrix_t, (K - 1)> C_bar;
    array<Model::state_vector_t, (K - 1)> Sigma_bar;
    array<Model::state_vector_t, (K - 1)> z_bar;

    print("Initializing solver.\n");
    optimization_problem::SecondOrderConeProgram socp = sc::build_successive_convexification_SOCP(model,
                                                                                                  weight_trust_region_sigma, weight_trust_region_xu, weight_virtual_control,
                                                                                                  X, U, sigma, A_bar, B_bar, C_bar, Sigma_bar, z_bar);

    // Cache indices for performance
    const size_t sigma_index = socp.get_tensor_variable_index("sigma", {});
    size_t X_indices[Model::state_dim_][K];
    size_t U_indices[Model::input_dim_][K];
    for (size_t k = 0; k < K; k++)
    {
        for (size_t i = 0; i < Model::state_dim_; ++i)
        {
            X_indices[i][k] = socp.get_tensor_variable_index("X", {i, k});
        }
        for (size_t i = 0; i < Model::input_dim_; ++i)
        {
            U_indices[i][k] = socp.get_tensor_variable_index("U", {i, k});
        }
    }

    EcosWrapper solver(socp);

    print("Starting Successive Convexification.\n");
    const size_t iterations = 50;
    const double timer_total = tic();
    for (size_t it = 0; it < iterations; it++)
    {
        string itString = format("<Iteration {}>", it);
        print("{:=^{}}\n", itString, 60);

        weight_trust_region_xu *= 2.;

        const double timer_iteration = tic();
        double timer = tic();
        calculate_discretization(model, sigma, X, U, A_bar, B_bar, C_bar, Sigma_bar, z_bar);
        print("{:<{}}{:.2f}ms\n", "Time, discretization:", 50, toc(timer));

        timer = tic();
        string file_name_prefix = format("{}iteration{:03}_", get_output_path(), it);

        // Write solution to files
        {
            ofstream f(file_name_prefix + "X.txt");
            f << X;
        }
        {
            ofstream f(file_name_prefix + "U.txt");
            f << U;
        }
        print("{:<{}}{:.2f}ms\n", "Time, solution file:", 50, toc(timer));

        // Write problem to file
        // timer = tic();
        // {
        //     ofstream f(file_name_prefix + "problem.txt");
        //     socp.print_problem(f);
        // }
        // print("{:<{}}{:.2f}ms\n", "Time, problem file:", 50, toc(timer));

        // save matrices for debugging
        // if (it == 0)
        // {
        //     for (unsigned int k = 0; k < K - 1; k++)
        //     {
        //         {
        //             ofstream f(format("{}z_bar{}.txt", get_output_path(), k));
        //             f << z_bar.at(k);
        //         }
        //     }
        // }

        timer = tic();
        solver.solve_problem();
        print("{:<{}}{:.2f}ms\n", "Time, solver:", 50, toc(timer));

        timer = tic();
        if (!socp.feasibility_check(solver.get_solution_vector()))
        {
            print("ERROR: Solver produced an invalid solution.\n");
            return EXIT_FAILURE;
        }
        print("{:<{}}{:.2f}ms\n", "Time, solution check:", 50, toc(timer));

        // Read solution
        for (size_t k = 0; k < K; k++)
        {
            for (size_t i = 0; i < Model::state_dim_; ++i)
            {
                X(i, k) = solver.get_solution_value(X_indices[i][k]);
            }
            for (size_t i = 0; i < Model::input_dim_; ++i)
            {
                U(i, k) = solver.get_solution_value(U_indices[i][k]);
            }
        }
        sigma = solver.get_solution_value(sigma_index);

        print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
        print("\n");
        print("{:<{}}{: .4f}\n", "sigma", 50, sigma);
        print("{:<{}}{: .4f}\n", "norm1_nu", 50, solver.get_solution_value("norm1_nu", {}));
        print("{:<{}}{: .4f}\n", "Delta_sigma", 50, solver.get_solution_value("Delta_sigma", {}));
        print("{:<{}}{: .4f}\n", "norm2_Delta", 50, solver.get_solution_value("norm2_Delta", {}));
        print("\n");
        print("{:<{}}{:.2f}ms\n", "Time, iteration:", 50, toc(timer_iteration));
        print("\n");
        if (solver.get_solution_value("norm2_Delta", {}) < delta_tol && solver.get_solution_value("norm1_nu", {}) < nu_tol)
        {
            print("Converged after {} iterations.\n", it + 1);
            break;
        }
        else if (it == iterations - 1)
        {
            print("No convergence after {} iterations.\n", it + 1);
        }
    }
    print("{:<{}}{:.2f}ms\n", "Time, total:", 50, toc(timer_total));
}