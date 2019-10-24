#pragma once

#include "activeModel.hpp"

namespace scpp::discretization
{

struct Data
{
    Model::state_matrix_v_t A;
    Model::control_matrix_v_t B;
    Model::control_matrix_v_t C;
    Model::state_vector_v_t s;
    Model::state_vector_v_t z;

    void initialize(size_t K, bool interpolate_input, bool free_final_time);

    bool interpolatedInput() const;

    bool variableTime() const;
};

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
    Data &dd);

} // namespace scpp::discretization