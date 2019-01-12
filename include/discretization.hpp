#pragma once

#include <array>
#include <boost/numeric/odeint.hpp>

#include "eigenIntegration.hpp"
#include "active_model.hpp"

using std::array;

void calculate_discretization(
    Model &model,
    double &sigma,
    Eigen::Matrix<double, Model::state_dim_, K> &X,
    Eigen::Matrix<double, Model::input_dim_, K> &U,
    array<Model::state_matrix_t, (K - 1)> &A_bar,
    array<Model::state_input_matrix_t, (K - 1)> &B_bar,
    array<Model::state_input_matrix_t, (K - 1)> &C_bar,
    array<Model::state_vector_t, (K - 1)> &Sigma_bar,
    array<Model::state_vector_t, (K - 1)> &z_bar);