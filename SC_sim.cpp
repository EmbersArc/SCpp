#include <experimental/filesystem>

#include "SCAlgorithm.hpp"
#include "timing.hpp"

namespace fs = std::experimental::filesystem;

fs::path getOutputPath() { return fs::path("..") / "output" / Model::getModelName(); }

int main()
{
    auto model = std::make_shared<Model>();
    model->loadParameters();

    SCAlgorithm solver(model);

    solver.initialize();

    const size_t sim_steps = 1;

    Model::state_vector_v_t X;
    Model::input_vector_v_t U;
    Model::state_vector_v_t X_sim(sim_steps);
    Model::input_vector_v_t U_sim(sim_steps);

    double t;

    double timer_run = tic();
    for (size_t i = 0; i < sim_steps; i++)
    {
        fmt::print("\nSIMULATION STEP {}:\n", i);
        solver.solve(i > 0);
        solver.getSolution(X, U, t);

        X_sim.push_back(X.at(0));
        U_sim.push_back(U.at(0));

        model->p.x_init = X.at(1);
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
        for (size_t i = 0; i < X_sim.size(); i++)
        {
            f << X_sim.at(i).transpose().format(CSVFormat) << "\n";
        }
    }
    {
        std::ofstream f(outputPath / "U_sim.txt");
        for (size_t i = 0; i < X_sim.size(); i++)
        {
            f << U_sim.at(i).transpose().format(CSVFormat) << "\n";
        }
    }
    {
        std::ofstream f(outputPath / "t_sim.txt");
        f << t / (X.size() - 1) * sim_steps;
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
}