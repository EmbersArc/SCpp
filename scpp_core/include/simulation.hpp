#pragma once

#include "activeModel.hpp"

namespace scpp
{

class ODE
{
public:
    ODE(Model::ptr_t model, double dt, const Model::input_vector_t &u0, const Model::input_vector_t &u1);
    void operator()(const Model::state_vector_t &f, Model::state_vector_t &dfdt, const double t);

private:
    Model::ptr_t model;
    Model::input_vector_t u0, u1;
    double dt;
};

void simulate(Model::ptr_t model, double dt,
              const Model::input_vector_t &u0,
              const Model::input_vector_t &u1,
              Model::state_vector_t &x);

} // namespace scpp