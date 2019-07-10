#include <experimental/filesystem>

#include "positionController.hpp"
#include "simulation.hpp"
#include "timing.hpp"

namespace fs = std::experimental::filesystem;

fs::path getOutputPath() { return fs::path("..") / "output" / Model::getModelName(); }

int main()
{
    auto model = std::make_shared<Model>();
    model->p.loadFromFile();
    model->initializeModel();

    const double sim_time = 15.;

    Model::dynamic_matrix_t X;
    Model::dynamic_matrix_t U;

    Model::state_vector_v_t X_sim;
    Model::input_vector_v_t U_sim;
    std::vector<double> times;

    Model::input_vector_t u(0, 0, -model->p.g_I.z() * model->p.m);
    Model::state_vector_t x = model->p.x_init;

    double current_time = 0.;
    size_t sim_step = 0;

    const double min_timestep = 0.010;

    const double t_start_run = tic();
    while (current_time < sim_time)
    {
        fmt::print("{:=^{}}\n", fmt::format("<SIMULATION STEP {}>", sim_step), 60);

        X_sim.push_back(x);
        U_sim.push_back(u);
        times.push_back(current_time);

        // move solve_time forward
        sim::simulate(*model, min_timestep, x, u, u, x);
        current_time += min_timestep;
    
        // get the calculated input
        u = woodPositionController(model, x, model->p.x_final);
        // u = customPositionController(model, x, model->p.x_final);

        sim_step++;

        // if ((x - model->p.x_final).head<3>().norm() < 0.01 and x.segment<3>(3).norm() < 0.01)
        // {
        //     break;
        // }
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