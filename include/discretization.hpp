#pragma once

#include <boost/numeric/odeint.hpp>

#include "eigenIntegration.hpp"
#include "activeModel.hpp"

void calculateDiscretization(
    Model &model,
    double &sigma,
    const Eigen::MatrixXd &X,
    const Eigen::MatrixXd &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &Sigma_bar,
    Model::state_vector_v_t &z_bar);