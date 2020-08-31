#include "sc_dynamic.hpp"

trajectory_data_t sc_dynamic(std::shared_ptr<Model> model)
{
    scpp::SCAlgorithm solver(model);

    solver.initialize();

    trajectory_data_t td;

    solver.solve();
    solver.getSolution(td);

    return td;
}