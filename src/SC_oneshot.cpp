#include <experimental/filesystem>

#include "SCAlgorithm.hpp"
#include "timing.hpp"

namespace fs = std::experimental::filesystem;

fs::path getOutputPath() { return fs::path("..") / "output" / Model::getModelName(); }

/**
 * @brief Computes a single SC trajectory and saves all intermediate iterations.
 * 
 */
int main()
{
    auto model = std::make_shared<Model>();
    model->loadParameters();

    SCAlgorithm solver(model);

    solver.initialize();

    std::vector<Model::state_vector_v_t> all_X;
    std::vector<Model::input_vector_v_t> all_U;
    std::vector<double> all_times;

    solver.solve();
    solver.getAllSolutions(all_X, all_U, all_times);

    // write solution to files
    double timer = tic();
    const fs::path outputPath = getOutputPath() / "SC" / std::to_string(size_t(timer));

    const Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                    Eigen::DontAlignCols,
                                    ", ", "\n");

    for (size_t k = 0; k < all_X.size(); k++)
    {
        const fs::path iterationPath = outputPath / std::to_string(k);
        if (not fs::exists(iterationPath) and not fs::create_directories(iterationPath))
        {
            throw std::runtime_error("Could not create output directory!");
        }

        {
            std::ofstream f(iterationPath / "X.txt");
            for (auto &state : all_X.at(k))
            {
                f << state.transpose().format(CSVFormat) << "\n";
            }
        }
        {
            std::ofstream f(iterationPath / "U.txt");
            for (auto &input : all_U.at(k))
            {
                f << input.transpose().format(CSVFormat) << "\n";
            }
        }
        {
            std::ofstream f(iterationPath / "t.txt");
            f << all_times.at(k);
        }
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
}