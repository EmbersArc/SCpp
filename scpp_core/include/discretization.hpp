#pragma once

#include "activeModel.hpp"
#include "discretizationData.hpp"

namespace scpp::discretization
{

void exactLinearDiscretization(Model::ptr_t model,
                               double ts,
                               const Model::state_vector_t &x_eq,
                               const Model::input_vector_t &u_eq,
                               Model::state_matrix_t &A,
                               Model::control_matrix_t &B,
                               Model::state_vector_t &z);

void multipleShooting(
    Model::ptr_t model,
    double time,
    const Model::state_vector_v_t &X,
    const Model::input_vector_v_t &U,
    DiscretizationData &dd);

} // namespace scpp::discretization