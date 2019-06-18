#pragma once

#include "activeModelMPC.hpp"

void eulerLinearDiscretization(Model &model,
                               double &ts,
                               const Model::state_vector_t &x_eq,
                               const Model::input_vector_t &u_eq,
                               Model::state_matrix_t &A,
                               Model::control_matrix_t &B);

void exactLinearDiscretization(Model &model,
                               double &ts,
                               const Model::state_vector_t &x_eq,
                               const Model::input_vector_t &u_eq,
                               Model::state_matrix_t &A,
                               Model::control_matrix_t &B);
