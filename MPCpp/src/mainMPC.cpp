#include <filesystem>

#include "MPCAlgorithm.hpp"
#include "timing.hpp"

namespace fs = std::filesystem;

fs::path getOutputPath()
{
    return fs::path("..") / "output" / Model::getModelName();
}

int main()
{
    auto model = std::make_shared<Model>();
    MPCAlgorithm solver(model);

    double T;
    model->getTimeHorizon(T);
    const size_t sim_steps = 100;

    solver.initialize();
    Model::dynamic_matrix_t X;
    Model::dynamic_matrix_t U;
    Model::dynamic_matrix_t X_sim(size_t(Model::state_dim), sim_steps);
    Model::dynamic_matrix_t U_sim(size_t(Model::input_dim), sim_steps);

    solver.setDesiredState(model->p.x_final);
    solver.setInitialState(model->p.x_init);

    double timer_run = tic();
    for (size_t i = 0; i < sim_steps; i++)
    {
        fmt::print("\nSIMULATION STEP {}:\n", i);
        solver.solve();
        solver.getSolution(X, U);

        X_sim.col(i) = X.col(0);
        U_sim.col(i) = U.col(0);

        model->p.x_init = X.col(1);
        solver.setInitialState(model->p.x_init);
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

    {
        std::ofstream f(outputPath / "X_sim.txt");
        f << X_sim;
    }
    {
        std::ofstream f(outputPath / "U_sim.txt");
        f << U_sim;
    }
    {
        std::ofstream f(outputPath / "t_sim.txt");
        f << T / (X.cols() - 1) * sim_steps;
    }
    {
        std::ofstream f(outputPath / "X.txt");
        f << X;
    }
    {
        std::ofstream f(outputPath / "U.txt");
        f << U;
    }
    {
        std::ofstream f(outputPath / "t.txt");
        f << T;
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
}