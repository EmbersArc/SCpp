#include <experimental/filesystem>

#include "LQRAlgorithm.hpp"
#include "simulation.hpp"
#include "timing.hpp"
#include "vectorOperations.hpp"

namespace fs = std::experimental::filesystem;

fs::path getOutputPath() { return fs::path("..") / "output" / Model::getModelName(); }

int main()
{
    auto model = std::make_shared<Model>();
    model->loadParameters();
    model->initializeModel();

    lqr::LQRAlgorithm solver(model);

    const double sim_time = 15.;
    const double write_steps = 30;

    solver.initialize();

    Model::state_vector_v_t X;
    Model::input_vector_v_t U;

    Model::state_vector_v_t X_sim;
    Model::input_vector_v_t U_sim;
    std::vector<double> t_sim;

    Model::input_vector_t u;
    Model::state_vector_t x = model->p.x_init;

    solver.setFinalState(model->p.x_final);
    std::cout << model->p.x_init << "\n";
    std::cout << model->p.x_final << "\n";

    double t = 0.;
    size_t sim_step = 0;

    const double time_step = 0.010;

    const double run_timer = tic();
    while (t < sim_time)
    {
        fmt::print("{:=^{}}\n", fmt::format("<SIMULATION STEP {}>", sim_step), 60);

        solver.setInitialState(x);

        // solve with current state
        solver.solve();

        // get the calculated input
        solver.getSolution(u);

        u.z() = std::max(model->p.T_min, u.z());

        // constrain input
        const double c = std::tan(model->p.gimbal_max) * u.z();
        if (u.head<2>().norm() > c)
        {
            u.head<2>() = c * u.head<2>().normalized();
        }
        if (u.norm() > model->p.T_max)
        {
            u = model->p.T_max * u.normalized();
        }

        // move time forward
        sim::simulate(model, time_step, x, u, u, x);
        t += time_step;

        X_sim.push_back(x);
        U_sim.push_back(u);
        t_sim.push_back(t);

        sim_step++;

        if ((x - model->p.x_final).norm() < 0.02)
        {
            break;
        }
    }
    fmt::print("\n");
    fmt::print("{:=^{}}\n", fmt::format("<SIMULATION FINISHED>"), 60);
    fmt::print("{:<{}}{:.2f}s\n", "Runtime:", 50, 0.001 * toc(run_timer));
    fmt::print("{:<{}}{:.2f}s\n", "Simulated time:", 50, t);
    const double freq = double(sim_step) / t;
    fmt::print("{:<{}}{:.2f}Hz\n", "Average frequency:", 50, freq);
    fmt::print("\n");

    // write solution to files
    double write_timer = tic();
    fs::path outputPath = getOutputPath() / "LQR" / std::to_string(size_t(write_timer)) / "0";
    if (not fs::exists(outputPath) and not fs::create_directories(outputPath))
    {
        throw std::runtime_error("Could not create output directory!");
    }

    const Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                    Eigen::DontAlignCols,
                                    ", ", "\n");

    X_sim = reduce_vector(X_sim, write_steps);
    U_sim = reduce_vector(U_sim, write_steps);
    t_sim = reduce_vector(t_sim, write_steps);
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