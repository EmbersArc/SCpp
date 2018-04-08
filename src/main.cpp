/*
 * 
 * Independent implementation of 
 *     Successive Convexification for 6-DoF Mars Rocket Powered Landing 
 *     with Free-Final-Time (Michael Szmuk, Behcet Acikmese)
 * 
 * https://arxiv.org/abs/1802.03827
 * 
 */

#include "active_model.hpp"
#include "EcosWrapper.hpp"
#include "MosekWrapper.hpp"
#include "Discretization.hpp"
#include "SuccessiveConvexificationSOCP.hpp"

#include <iostream>
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using std::array;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::setw;
using std::setfill;



string get_output_path() {
    return "../output/" + Model::get_name() + "/";
}

void make_output_path() {
    string path = "mkdir -p " + get_output_path();
    system(path.c_str());
}


int main() {
    make_output_path();
    Model model;

    // trajectory points
    const double dt = 1 / double(K-1);

    double weight_trust_region_sigma = 1e-1;
    double weight_trust_region_xu = 1e-1;
    double weight_virtual_control = 10000;

    const size_t n_states = Model::n_states;
    const size_t n_inputs = Model::n_inputs;

    Eigen::Matrix<double, n_states, K> X;
    Eigen::Matrix<double, n_inputs, K> U;


    // START INITIALIZATION
    cout << "Starting initialization." << endl;
    model.initialize(X, U);
    cout << "Initialization finished." << endl;

    // START SUCCESSIVE CONVEXIFICATION
    
    double sigma = model.total_time_guess();


    array<Model::StateMatrix,   (K-1)> A_bar;
    array<Model::ControlMatrix, (K-1)> B_bar;
    array<Model::ControlMatrix, (K-1)> C_bar;
    array<Model::StateVector,   (K-1)> Sigma_bar;
    array<Model::StateVector,   (K-1)> z_bar;




    optimization_problem::SecondOrderConeProgram socp = build_successive_convexification_SOCP ( 
        model, weight_trust_region_sigma, weight_trust_region_xu, weight_virtual_control, X, U, sigma, A_bar, B_bar, C_bar, Sigma_bar, z_bar );


    // Cache indices for performance
    const size_t sigma_index = socp.get_tensor_variable_index("sigma", {});
    size_t X_indices[n_states][K];
    size_t U_indices[n_inputs][K];
    for (size_t k = 0; k < K; k++) {
        for (size_t i = 0; i < n_states; ++i) X_indices[i][k] = socp.get_tensor_variable_index("X",{i,k});
        for (size_t i = 0; i < n_inputs; ++i) U_indices[i][k] = socp.get_tensor_variable_index("U",{i,k});
    }


    //EcosWrapper solver(socp);
    MosekWrapper solver(socp);

    const size_t iterations = 40;
    for(size_t it = 0; it < iterations; it++) {
        clock_t begin_time = clock();

        weight_trust_region_xu *= 1.2;

        cout << "Iteration " << it << endl;
        cout << "Calculating new transition matrices." << endl;



        calculate_discretization ( model, sigma, X, U, A_bar, B_bar, C_bar, Sigma_bar, z_bar );


        // Write problem to file

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

        /************************************************************************************/
        cout << "Solving problem" << endl;
        begin_time = clock();
        solver.solve_problem();
        cout << endl << "Solver time: " << double( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds." << endl;
        if(!socp.feasibility_check(solver.get_solution_vector())) {
            cout << "ERROR: Solver produced an invalid solution." << endl;
            return EXIT_FAILURE;
        }

        // Read solution
        for (size_t k = 0; k < K; k++) {
            for (size_t i = 0; i < n_states; ++i) X(i,k) = solver.get_solution_value(X_indices[i][k]);
            for (size_t i = 0; i < n_inputs; ++i) U(i,k) = solver.get_solution_value(U_indices[i][k]);
        }
        sigma = solver.get_solution_value(sigma_index);


        // Write solution to files

        {
            ofstream f(file_name_prefix + "X.txt");
            f << X;
        }
        {
            ofstream f(file_name_prefix + "U.txt");
            f << U;
        }

        cout << "sigma   " << sigma << endl;
        cout << "norm2_nu   " << solver.get_solution_value("norm2_nu", {}) << endl;
        cout << "Delta_sigma   " << solver.get_solution_value("Delta_sigma", {}) << endl;
        cout << "norm2_Delta   " << solver.get_solution_value("norm2_Delta", {}) << endl;

        cout << "Iteration time: " << (double( clock () - begin_time ) /  CLOCKS_PER_SEC*1000.) << " msec" << endl;
        cout << "==========================================================" << endl;
    }
}