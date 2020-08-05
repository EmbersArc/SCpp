#include <experimental/filesystem>

#include "SCAlgorithm.hpp"
#include "timing.hpp"
#include "commonFunctions.hpp"

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

    scpp::SCAlgorithm solver(model);

    solver.initialize();

    std::vector<trajectory_data_t> all_td;

    solver.solve();
    solver.getAllSolutions(all_td);

    // write solution to files
    double timer = tic();
    const fs::path outputPath = getOutputPath() / "SC" / scpp::getTimeString();

    const Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                    Eigen::DontAlignCols,
                                    ", ", "\n");

    for (size_t k = 0; k < all_td.size(); k++)
    {
        const fs::path iterationPath = outputPath / std::to_string(k);
        if (not fs::exists(iterationPath) and not fs::create_directories(iterationPath))
        {
            throw std::runtime_error("Could not create output directory!");
        }

        {
            std::ofstream f(iterationPath / "X.txt");
            for (auto &state : all_td.at(k).X)
            {
                f << state.transpose().format(CSVFormat) << "\n";
            }
        }
        {
            std::ofstream f(iterationPath / "U.txt");
            for (auto &input : all_td.at(k).U)
            {
                f << input.transpose().format(CSVFormat) << "\n";
            }
        }
        {
            std::ofstream f(iterationPath / "t.txt");
            f << all_td.at(k).t;
        }
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
}