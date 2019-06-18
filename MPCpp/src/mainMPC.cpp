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

    const size_t sim_steps = 2;

    mpcSolver.initialize();
    Model::dynamic_matrix_t X;
    Model::dynamic_matrix_t X_full(size_t(Model::state_dim), sim_steps);
    Model::dynamic_matrix_t U;

    Model::state_vector_t x_init = model->par.x_init;

    mpcSolver.setDesiredState(model->par.x_final);

    for (size_t i = 0; i < sim_steps; i++)
    {
        X_full.col(i) = x_init;
        std::cout << x_init.transpose() << std::endl;
        mpcSolver.setInitialState(x_init);
        mpcSolver.solve();
        mpcSolver.getSolution(X, U);

        x_init = X.col(1);
    }

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
        std::ofstream f(outputPath / "X.txt");
        f << X;
    }
    {
        std::ofstream f(outputPath / "U.txt");
        f << U;
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution file:", 50, toc(timer));
}