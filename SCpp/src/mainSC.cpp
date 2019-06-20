#include <experimental/filesystem>

#include "SCAlgorithm.hpp"
#include "timing.hpp"

namespace fs = std::experimental::filesystem;

fs::path getOutputPath()
{
    return fs::path("..") / "output" / Model::getModelName();
}

int main()
{
    auto model = std::make_shared<Model>();
    SCAlgorithm scSolver(model);

    scSolver.initialize();
    scSolver.solve();

    Model::dynamic_matrix_t X;
    Model::dynamic_matrix_t U;
    double t;
    scSolver.getSolution(X, U, t);

    // write solution to files
    double timer = tic();
    fs::path outputPath = getOutputPath() / std::to_string(time(0));
    if (not fs::exists(outputPath) and not fs::create_directories(outputPath))
    {
        throw std::runtime_error("Could not create output directory!");
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
        f << t;
    }

    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
}