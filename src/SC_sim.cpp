#include <experimental/filesystem>

#include "SCAlgorithm.hpp"
#include "simulation.hpp"
#include "timing.hpp"
#include "vectorOperations.hpp"

namespace fs = std::experimental::filesystem;

fs::path getOutputPath() { return fs::path("..") / "output" / Model::getModelName(); }

Model::input_vector_t interpolatedInput(const Model::input_vector_v_t &U, double t, double total_time)
{
    const size_t K = U.size();
    const double time_step = total_time / (K - 1);
    const size_t i = std::min(size_t(t / time_step), K - 2);
    const Model::input_vector_t u0 = U.at(i);
    const Model::input_vector_t u1 = U.at(i + 1);

    const double t_intermediate = std::fmod(t, time_step) / time_step;
    const Model::input_vector_t u = u0 + (u1 - u0) * t_intermediate;
    return u;
}

/**
 * @brief Simulates a trajectory with the SC controller.
 * 
 * (Kind of broken right now)
 */
int main()
{
    auto model = std::make_shared<Model>();
    model->loadParameters();

    scpp::SCAlgorithm solver(model);

    solver.initialize();

    const double time_step = 0.015;

    Model::state_vector_v_t X;
    Model::input_vector_v_t U;
    Model::state_vector_v_t X_sim;
    Model::input_vector_v_t U_sim;

    double t;
    Model::state_vector_t &x = model->p.x_init;

    double timer_run = tic();
    size_t sim_step = 0;
    while ((x - model->p.x_final).norm() > 0.02 and sim_step < 100)
    {
        fmt::print("\nSIMULATION STEP {}:\n", sim_step);

        const bool warm_start = sim_step > 0;
        solver.solve(warm_start);
        solver.getSolution(X, U, t);

        const Model::input_vector_t u0 = U.at(0);
        const Model::input_vector_t u1 = interpolatedInput(U, time_step, t);

        scpp::simulate(model, time_step, x, u0, u1, x);

        X_sim.push_back(x);
        U_sim.push_back(u0);

        sim_step++;
    }
    fmt::print("\n");
    fmt::print("{:<{}}{:.2f}ms\n", fmt::format("Time, {} steps:", sim_step), 50, toc(timer_run));
    const double freq = double(sim_step) / (0.001 * toc(timer_run));
    fmt::print("{:<{}}{:.2f}Hz\n", "Average frequency:", 50, freq);
    fmt::print("\n");

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
        f << t / (X.size() - 1) * sim_step;
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
}