#include "activeModel.hpp"

#include "Eigen/Dense"

bool ComputeLQR(const Model::state_matrix_t &Q,
                const Model::input_matrix_t &R,
                const Model::state_matrix_t &A,
                const Model::control_matrix_t &B,
                Model::feedback_matrix_t &K,
                bool RisDiagonal);