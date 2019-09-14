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

    Model::state_vector_v_t X;
    Model::input_vector_v_t U;
    double t;

    double timer_run = tic();

    solver.solve();
    solver.getSolution(X, U, t);

    fmt::print("\n");
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution:", 50, toc(timer_run));
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
        std::ofstream f(outputPath / "X.txt");
        for (size_t i = 0; i < X.size(); i++)
        {
            f << X.at(i).transpose().format(CSVFormat) << "\n";
        }
    }
    {
        std::ofstream f(outputPath / "U.txt");
        for (size_t i = 0; i < U.size(); i++)
        {
            f << U.at(i).transpose().format(CSVFormat) << "\n";
        }
    }
    {
        std::ofstream f(outputPath / "t.txt");
        f << t;
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
}