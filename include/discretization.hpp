#pragma once

#include <array>
#include <boost/numeric/odeint.hpp>

#include "eigenIntegration.hpp"
#include "active_model.hpp"

using std::array;

void calculate_discretization(
    Model &model,
    double &sigma,
    Eigen::MatrixXd &X,
    Eigen::MatrixXd &U,
    vector<Model::state_matrix_t> &A_bar,
    vector<Model::control_matrix_t> &B_bar,
    vector<Model::control_matrix_t> &C_bar,
    vector<Model::state_vector_t> &Sigma_bar,
    vector<Model::state_vector_t> &z_bar);