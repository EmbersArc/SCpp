/*
 * 
 * Independent implementation of 
 *     Successive Convexification for 6-DoF Mars Rocket Powered Landing 
 *     with Free-Final-Time (Michael Szmuk, Behcet Acikmese)
 * 
 * https://arxiv.org/abs/1802.03827
 * 
 */

#include <iostream>
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "active_model.hpp"
#include "ecosWrapper.hpp"
#include "mosekWrapper.hpp"
#include "discretization.hpp"
#include "successiveConvexificationSOCP.hpp"
#include "timing.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using std::array;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::setfill;
using std::setw;

string get_output_path()
{
    return "../output/" + Model::get_name() + "/";
}

void clear_output_path()
{
    string command = "rm -r " + get_output_path();
    assert(system(command.c_str()) == 0);
}

void make_output_path()
{
    string command = "mkdir -p " + get_output_path();
    assert(system(command.c_str()) == 0);
}

int main()
{
    clear_output_path();
    make_output_path();
    Model model;

    double weight_trust_region_sigma = 1e-1;
    double weight_trust_region_xu = 1e-3;
    double weight_virtual_control = 1e5;

    double nu_tol = 1e-8;
    double delta_tol = 1e-3;

    const size_t n_states = Model::n_states;
    const size_t n_inputs = Model::n_inputs;

    Eigen::Matrix<double, n_states, K> X;
    Eigen::Matrix<double, n_inputs, K> U;

    model.initialize(X, U);

    double sigma = model.total_time_guess();

    array<Model::StateMatrix, (K - 1)> A_bar;
    array<Model::ControlMatrix, (K - 1)> B_bar;
    array<Model::ControlMatrix, (K - 1)> C_bar;
    array<Model::StateVector, (K - 1)> Sigma_bar;
    array<Model::StateVector, (K - 1)> z_bar;

    optimization_problem::SecondOrderConeProgram socp = build_successive_convexification_SOCP(
        model, weight_trust_region_sigma, weight_trust_region_xu, weight_virtual_control, X, U, sigma, A_bar, B_bar, C_bar, Sigma_bar, z_bar);

    // Cache indices for performance
    const size_t sigma_index = socp.get_tensor_variable_index("sigma", {});
    size_t X_indices[n_states][K];
    size_t U_indices[n_inputs][K];
    for (size_t k = 0; k < K; k++)
    {
        for (size_t i = 0; i < n_states; ++i)
            X_indices[i][k] = socp.get_tensor_variable_index("X", {i, k});
        for (size_t i = 0; i < n_inputs; ++i)
            U_indices[i][k] = socp.get_tensor_variable_index("U", {i, k});
    }

    EcosWrapper solver(socp);
    //    MosekWrapper solver(socp);

    const size_t iterations = 30;

    const double timer_total = tic();
    for (size_t it = 0; it < iterations; it++)
    {

        //        weight_trust_region_xu *= 1.2;

        const double timer_iteration = tic();
        double timer = tic();
        calculate_discretization(model, sigma, X, U, A_bar, B_bar, C_bar, Sigma_bar, z_bar);
        cout << "Time, discretization: " << toc(timer) << " ms" << endl;

        //        // Write problem to file
        // timer = tic();
        string file_name_prefix;
        {
            ostringstream file_name_prefix_ss;
            file_name_prefix_ss << get_output_path() << "iteration"
                                << setfill('0') << setw(3) << it << "_";
            file_name_prefix = file_name_prefix_ss.str();
        }

        {
            ofstream f(file_name_prefix + "problem.txt");
            socp.print_problem(f);
        }
        cout << "Time, problem file: " << toc(timer) << " ms" << endl;

        // save matrices for debugging
        if (it == 0)
        {
            for (unsigned int k = 0; k < K - 1; k++)
            {
                {
                    ofstream f(get_output_path() + "z_bar" + std::to_string(k) + ".txt");
                    f << z_bar.at(k);
                }
            }
        }

        timer = tic();
        solver.solve_problem();
        cout << "Time, solver: " << toc(timer) << " ms" << endl;

        timer = tic();
        if (!socp.feasibility_check(solver.get_solution_vector()))
        {
            cout << "ERROR: Solver produced an invalid solution." << endl;
            return EXIT_FAILURE;
        }
        cout << "Time, solution check: " << toc(timer) << " ms" << endl;

        // Read solution
        for (size_t k = 0; k < K; k++)
        {
            for (size_t i = 0; i < n_states; ++i)
                X(i, k) = solver.get_solution_value(X_indices[i][k]);
            for (size_t i = 0; i < n_inputs; ++i)
                U(i, k) = solver.get_solution_value(U_indices[i][k]);
        }
        sigma = solver.get_solution_value(sigma_index);

        //         Write solution to files
        timer = tic();
        {
            ofstream f(file_name_prefix + "X.txt");
            f << X;
        }
        {
            ofstream f(file_name_prefix + "U.txt");
            f << U;
        }
        cout << "Time, solution files: " << toc(timer) << " ms" << endl;

        cout << "sigma   " << sigma << endl;
        cout << "norm1_nu   " << solver.get_solution_value("norm1_nu", {}) << endl;
        cout << "Delta_sigma   " << solver.get_solution_value("Delta_sigma", {}) << endl;
        cout << "norm2_Delta   " << solver.get_solution_value("norm2_Delta", {}) << endl;
        cout << "Time, iteration: " << toc(timer_iteration) << " ms" << endl;
        cout << "==========================================================" << endl;

        if (solver.get_solution_value("norm2_Delta", {}) < delta_tol && solver.get_solution_value("norm1_nu", {}) < nu_tol)
        {
            cout << "Converged after " << it << " iterations." << endl;
            break;
        }
        else if (it == iterations - 1)
        {
            cout << "No convergence after " << iterations << " iterations." << endl;
        }
    }
    cout << "Time, total: " << toc(timer_total) << " ms" << endl;
}