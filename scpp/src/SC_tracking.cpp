#include <experimental/filesystem>

#include "SCAlgorithm.hpp"
#include "LQRTracker.hpp"
#include "timing.hpp"
#include "simulation.hpp"
#include "commonFunctions.hpp"

namespace fs = std::experimental::filesystem;

fs::path getOutputPath() { return fs::path("..") / "output" / Model::getModelName(); }

/**
 * @brief Computes a single SC trajectory and tracks it via LQR.
 * 
 */
int main()
{
    auto model = std::make_shared<Model>();
    model->loadParameters();

    scpp::SCAlgorithm solver(model);

    solver.initialize();

    trajectory_data_t td;

    solver.solve();
    solver.getSolution(td);

    // calculate LQR gains
    fmt::print("\n");
    const double gains_timer = tic();
    scpp::LQRTracker tracker(model, td);
    fmt::print("{:<{}}{:.2f}ms\n", "Time, LQR gains:", 50, 1 * toc(gains_timer));

    const double timestep = 0.01;

    // start simulation
    const double run_timer = tic();

    Model::state_vector_v_t X_sim;
    Model::input_vector_v_t U_sim;
    std::vector<double> t_sim;
    const double t_max = td.t;
    const double write_steps = 30;

    Model::state_vector_t x = model->p.x_init;

    double t = 0.;
    size_t sim_step = 0;

    while (t < t_max)
    {
        // get the calculated input
        Model::input_vector_t u;
        tracker.getInput(t, x, u);

        // move solve_time forward
        scpp::simulate(model, timestep, u, u, x);
        t += timestep;

        X_sim.push_back(x);
        U_sim.push_back(u);
        t_sim.push_back(t);

        sim_step++;

        if (x.hasNaN())
        {
            throw std::runtime_error("State has NaN.");
        }

        if ((x - model->p.x_final).norm() < 0.01)
        {
            break;
        }
    }

    fmt::print("\n");
    fmt::print("Simulating trajectory.\n");
    fmt::print("Finished after {} steps.\n", sim_step + 1);
    fmt::print("Final error: {:.2f}.\n", (x - model->p.x_final).norm());
    fmt::print("{:<{}}{:.2f}ms\n", "Time, simulation:", 50, toc(run_timer));
    fmt::print("{:<{}}{:.2f}s\n", "Simulated time:", 50, t);
    fmt::print("\n");

    // write solution to files
    double write_timer = tic();
    fs::path outputPath = getOutputPath() / "SC_tracking" / scpp::getTimeString() / "0";
    if (not fs::exists(outputPath) and not fs::create_directories(outputPath))
    {
        throw std::runtime_error("Could not create output directory!");
    }

    const Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                    Eigen::DontAlignCols,
                                    ", ", "\n");

    X_sim = scpp::reduce_vector(X_sim, write_steps);
    U_sim = scpp::reduce_vector(U_sim, write_steps);
    t_sim = scpp::reduce_vector(t_sim, write_steps);
    {
        std::ofstream f(outputPath / "X.txt");
        for (auto &state : X_sim)
        {
            f << state.transpose().format(CSVFormat) << "\n";
        }
    }
    {
        std::ofstream f(outputPath / "U.txt");
        for (auto &input : U_sim)
        {
            f << input.transpose().format(CSVFormat) << "\n";
        }
    }
    {
        std::ofstream f(outputPath / "t.txt");
        for (auto &time : t_sim)
        {
            f << time << "\n";
        }
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(write_timer));
}