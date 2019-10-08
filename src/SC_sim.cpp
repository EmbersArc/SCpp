#include <experimental/filesystem>

#include "SCAlgorithm.hpp"
#include "simulation.hpp"
#include "timing.hpp"

namespace fs = std::experimental::filesystem;

fs::path getOutputPath() { return fs::path("..") / "output" / Model::getModelName(); }

/**
 * @brief Simulates a trajectory with the SC controller.
 * 
 * (Kind of broken right now)
 */
int main()
{
    auto model = std::make_shared<Model>();
    model->loadParameters();

    SCAlgorithm solver(model);

    solver.initialize();

    const size_t sim_steps = 100;
    const double time_step = 0.01;

    Model::state_vector_v_t X;
    Model::input_vector_v_t U;
    Model::state_vector_v_t X_sim(sim_steps);
    Model::input_vector_v_t U_sim(sim_steps);

    double t;
    Model::state_vector_t x = model->p.x_init;

    double timer_run = tic();
    for (size_t i = 0; i < sim_steps; i++)
    {
        fmt::print("\nSIMULATION STEP {}:\n", i);
        solver.solve(i > 0);
        solver.getSolution(X, U, t);

        // move solve_time forward
        const Model::input_vector_t u0 = U.at(0);
        const Model::input_vector_t u1 = U.at(1);

        sim::simulate(model, time_step, x, u0, u1, x);

        X_sim.push_back(x);
        U_sim.push_back(u0);

        if ((x - model->p.x_final).norm() < 0.02)
        {
            break;
        }
    }
    fmt::print("\n");
    fmt::print("{:<{}}{:.2f}ms\n", fmt::format("Time, {} steps:", sim_steps), 50, toc(timer_run));
    const double freq = double(sim_steps) / (0.001 * toc(timer_run));
    fmt::print("{:<{}}{:.2f}Hz\n", "Average frequency:", 50, freq);
    fmt::print("\n");

    // write solution to files
    double timer = tic();
    fs::path outputPath = getOutputPath() / std::to_string(0);
    if (not fs::exists(outputPath) and not fs::create_directories(outputPath))
    {
        throw std::runtime_error("Could not create output directory!");
    }

    const Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                    Eigen::DontAlignCols,
                                    ", ", "\n");

    {
        std::ofstream f(outputPath / "X_sim.txt");
        for (auto &x : X_sim)
        {
            f << x.transpose().format(CSVFormat) << "\n";
        }
    }
    {
        std::ofstream f(outputPath / "U_sim.txt");
        for (auto &u : U_sim)
        {
            f << u.transpose().format(CSVFormat) << "\n";
        }
    }
    {
        std::ofstream f(outputPath / "t_sim.txt");
        f << t / (X.size() - 1) * sim_steps;
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
}