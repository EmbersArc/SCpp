#pragma once

#include "activeModel.hpp"

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
    double T,
    const Model::state_vector_v_t &X,
    const Model::input_vector_v_t &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::state_vector_v_t &z_bar);

void multipleShooting(
    Model::ptr_t model,
    double T,
    const Model::state_vector_v_t &X,
    const Model::input_vector_v_t &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &z_bar);

void multipleShooting(
    Model::ptr_t model,
    double T,
    const Model::state_vector_v_t &X,
    const Model::input_vector_v_t &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::state_vector_v_t &S_bar,
    Model::state_vector_v_t &z_bar);

void multipleShooting(
    Model::ptr_t model,
    double T,
    const Model::state_vector_v_t &X,
    const Model::input_vector_v_t &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &S_bar,
    Model::state_vector_v_t &z_bar);

} // namespace scpp::discretization