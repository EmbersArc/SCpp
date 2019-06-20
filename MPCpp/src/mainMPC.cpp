#include <experimental/filesystem>

#include "MPCAlgorithm.hpp"
#include "timing.hpp"

namespace fs = std::experimental::filesystem;

fs::path getOutputPath()
{
    return fs::path("..") / "output" / Model::getModelName();
}

int main()
{
    auto model = std::make_shared<Model>();
    MPCAlgorithm mpcSolver(model);

    const size_t sim_steps = 200;

    mpcSolver.initialize();
    Model::dynamic_matrix_t X;
    Model::dynamic_matrix_t X_full(size_t(Model::state_dim), sim_steps);
    Model::dynamic_matrix_t U_full(size_t(Model::input_dim), sim_steps);
    Model::dynamic_matrix_t U;

    Model::state_vector_t x_init = model->par.x_init;

    mpcSolver.setDesiredState(model->par.x_final);

    double timer_run = tic();
    for (size_t i = 0; i < sim_steps; i++)
    {
        fmt::print("Step {}:\n", i);
        X_full.col(i) = x_init;
        mpcSolver.setInitialState(x_init);
        mpcSolver.solve();
        mpcSolver.getSolution(X, U);

        x_init = X.col(1);
        U_full.col(i) = U.col(0);
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
        std::ofstream f(outputPath / "X_full.txt");
        f << X_full;
    }
    {
        std::ofstream f(outputPath / "U_full.txt");
        f << U_full;
    }
    {
        std::ofstream f(outputPath / "X.txt");
        f << X;
    }
    {
        std::ofstream f(outputPath / "U.txt");
        f << U;
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
}