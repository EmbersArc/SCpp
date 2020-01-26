#pragma once

#include "activeModel.hpp"

namespace scpp
{

void simulate(Model::ptr_t model, double dt,
              const Model::input_vector_t &u0,
              const Model::input_vector_t &u1,
              Model::state_vector_t &x);

} // namespace scpp