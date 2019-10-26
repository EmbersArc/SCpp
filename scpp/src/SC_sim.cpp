#include <experimental/filesystem>

#include "SCAlgorithm.hpp"
#include "simulation.hpp"
#include "timing.hpp"
#include "vectorOperations.hpp"
#include "commonFunctions.hpp"

using fmt::format;
using fmt::print;
namespace fs = std::experimental::filesystem;

fs::path getOutputPath() { return fs::path("..") / "output" / Model::getModelName(); }

/**
 * @brief Simulates a trajectory with the SC controller.
 * 
 *
 */
int main()
{
    auto model = std::make_shared<Model>();
    model->loadParameters();

    scpp::SCvxAlgorithm solver(model);

    solver.initialize();

    const double time_step = 0.05;
    const size_t max_steps = 2;

    TrajectoryData td;

    Model::state_vector_v_t X_sim;
    Model::input_vector_v_t U_sim;

    Model::state_vector_t &x = model->p.x_init;

    double timer_run = tic();
    size_t sim_step = 0;
    while ((x - model->p.x_final).norm() > 0.02 and sim_step < max_steps)
    {
        print("\n{:*^{}}\n\n", format("<SIMULATION STEP {}>", sim_step), 60);

        const bool warm_start = sim_step > 0;
        solver.solve(warm_start);
        solver.getSolution(td);
        const size_t K = td.n_X();

        const Model::input_vector_t u0 = td.U.at(0);
        const Model::input_vector_t u1 = scpp::interpolatedInput(td.U, time_step, td.t, false);

        scpp::simulate(model, td.t / (K - 1), u0, u1, x);

        X_sim.push_back(x);
        U_sim.push_back(u0);

        sim_step++;
        // if (sim_step == max_steps)
        // {
        //     X_sim = reduce_vector(X_sim, X_sim.size() / 10);
        //     U_sim = reduce_vector(U_sim, U_sim.size() / 10);
        //     // add the planned trajectory
        //     X_sim.insert(X_sim.end(), X.begin(), X.end());
        //     U_sim.insert(U_sim.end(), U.begin(), U.end());
        // }
    }

    print("\n");
    print("{:<{}}{:.2f}ms\n", fmt::format("Time, {} steps:", sim_step), 50, toc(timer_run));
    const double freq = double(sim_step) / (0.001 * toc(timer_run));
    print("{:<{}}{:.2f}Hz\n", "Average frequency:", 50, freq);
    print("\n");

    // write solution to files
    double timer = tic();
    fs::path outputPath = getOutputPath() / "SC_sim" / std::to_string(size_t(timer)) / std::to_string(0);
    if (not fs::exists(outputPath) and not fs::create_directories(outputPath))
    {
        throw std::runtime_error("Could not create output directory!");
    }

    const Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                    Eigen::DontAlignCols,
                                    ", ", "\n");

    {
        std::ofstream f(outputPath / "X.txt");
        for (auto &state : X_sim)
        {
            f << state.transpose().format(CSVFormat) << "\n";
        }
    }
    {
        std::ofstream f(outputPath / "U.txt");
        for (auto &input : U_sim)
        {
            f << input.transpose().format(CSVFormat) << "\n";
        }
    }
    {
        std::ofstream f(outputPath / "t.txt");
        f << sim_step * time_step;
    }
    print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
}