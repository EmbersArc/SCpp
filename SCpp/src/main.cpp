#include "freeFinalTimeAlgorithm.hpp"
#include "timing.hpp"

namespace fs = std::filesystem;

std::string getOutputPath()
{
    return fmt::format("../output/{}", Model::getModelName());
}

int main()
{
    auto model = std::make_shared<Model>();
    freeFinalTimeAlgorithm scSolver(model);
    scSolver.initialize();
    scSolver.solve();

    Model::dynamic_matrix_t X;
    Model::dynamic_matrix_t U;
    double t;
    scSolver.getSolution(X, U, t);

    // write solution to files
    double timer = tic();
    std::string outputDirectory = fmt::format("{}/{}", getOutputPath(), time(0));

    if (not fs::exists(outputDirectory) and not fs::create_directories(outputDirectory))
    {
        throw std::runtime_error("Could not create output directory!");
    }

    {
        std::ofstream f(outputDirectory + "/X.txt");
        f << X;
    }
    {
        std::ofstream f(outputDirectory + "/U.txt");
        f << U;
    }
    {
        std::ofstream f(outputDirectory + "/t.txt");
        f << t;
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution file:", 50, toc(timer));
}