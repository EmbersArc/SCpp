#include <experimental/filesystem>

#include "SCAlgorithm.hpp"
#include "SCvxAlgorithm.hpp"
#include "timing.hpp"
#include "simulation.hpp"
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
    const fs::path outputPath = getOutputPath() / "SC" / std::to_string(size_t(timer));

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

    // verify by integrating
    // const auto &td = all_td.back();
    // Model::state_vector_t x = td.X.front();
    // const double timestep = td.t / (td.X.size() - 1);
    // const bool first_order_hold = td.interpolatedInput();

    // for (size_t k = 0; k < td.X.size() - 1; k++)
    // {
    //     scpp::simulate(model, timestep, td.U[k], first_order_hold ? td.U[k + 1] : td.U[k], x);
    // }

    // double full_error = (x - td.X.back()).lpNorm<1>();

    // fmt::print("\nVerifying solution...\n");
    // fmt::print("{:<{}}{:.5f}\n", "Final deviation:", 50, full_error);

    // save accelerations:
    // {
    //     const fs::path iterationPath = outputPath / "x";

    //     if (not fs::exists(iterationPath) and not fs::create_directories(iterationPath))
    //     {
    //         throw std::runtime_error("Could not create output directory!");
    //     }

    //     {
    //         std::ofstream f(iterationPath / "acc_passenger_b.txt");
    //         for (auto &pv : getAccelerationRotatingFrame(all_td.back(), Eigen::Vector3d(0., 0., 0.), 9.81))
    //         {
    //             f << pv.transpose().format(CSVFormat) << "\n";
    //         }
    //     }
    // }
}