#include <experimental/filesystem>

#include "MPCAlgorithm.hpp"
#include "simulation.hpp"
#include "timing.hpp"

namespace fs = std::experimental::filesystem;

fs::path getOutputPath() { return fs::path("..") / "output" / Model::getModelName(); }

// calculate an exponential moving average
double expMovingAverage(double previousAverage, double period, double newValue)
{
    const double factor = 2. / (period + 1.);
    const double result = (newValue - previousAverage) * factor + previousAverage;
    return result;
}

int main()
{
    auto model = std::make_shared<Model>();
    model->p.loadFromFile();

    mpc::MPCAlgorithm solver(model);

    const double sim_time = 15.;

    solver.initialize(true);

    Model::dynamic_matrix_t X;
    Model::dynamic_matrix_t U;

    Model::state_vector_v_t X_sim;
    Model::input_vector_v_t U_sim;
    std::vector<double> times;

    Model::input_vector_t u;
    u.setZero();
    Model::state_vector_t x = model->p.x_init;

    solver.setInitialState(x);
    solver.setFinalState(model->p.x_final);

    double current_time = 0.;
    size_t sim_step = 0;

    const double min_timestep = 0.010;
    double avg_solve_time = min_timestep;

    const double t_start_run = tic();
    while (current_time < sim_time)
    {
        fmt::print("{:=^{}}\n", fmt::format("<SIMULATION STEP {}>", sim_step), 60);

        // forward estimation
        Model::state_vector_t x_expected;
        sim::simulate(*model, avg_solve_time, x, u, u, x_expected);
        solver.setInitialState(x_expected);

        X_sim.push_back(x);
        U_sim.push_back(u);
        times.push_back(current_time);

        // solve with current state
        const double t_start_solve = tic();
        solver.solve();
        const double solve_time = std::max(toc(t_start_solve) * 0.001, min_timestep);
        avg_solve_time = expMovingAverage(avg_solve_time, 20, solve_time);
        fmt::print("{:<{}}{:.2f}ms\n", "Average solve time:", 50, avg_solve_time * 1000);

        // move solve_time forward
        sim::simulate(*model, solve_time, x, u, u, x);
        current_time += solve_time;

        // get the calculated input
        solver.getSolution(X, U);
        u = U.col(0);

        sim_step++;
    }
    fmt::print("\n");
    fmt::print("{:=^{}}\n", fmt::format("<SIMULATION FINISHED>"), 60);
    fmt::print("{:<{}}{:.2f}s\n", "Runtime:", 50, 0.001 * toc(t_start_run));
    fmt::print("{:<{}}{:.2f}s\n", "Simulation time:", 50, current_time);
    const double freq = double(sim_step) / current_time;
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
        std::ofstream f(outputPath / "X_sim.txt");
        Model::dynamic_matrix_t X;
        X.resize(Model::state_dim, sim_step);
        for (size_t i = 0; i < sim_step; i++)
        {
            X.col(i) = X_sim[i];
        }
        f << X;
    }
    {
        std::ofstream f(outputPath / "U_sim.txt");
        Model::dynamic_matrix_t U;
        U.resize(Model::input_dim, sim_step);
        for (size_t i = 0; i < sim_step; i++)
        {
            U.col(i) = U_sim[i];
        }
        f << U;
    }
    {
        std::ofstream f(outputPath / "t_sim.txt");
        for (size_t i = 0; i < sim_step; i++)
        {
            f << times[i] << '\n';
        }
    }
    fmt::print("{:<{}}{:.2f}ms\n", "Time, solution files:", 50, toc(timer));
}