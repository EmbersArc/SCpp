#pragma once

#include <array>
#include <boost/numeric/odeint.hpp>

#include "eigenIntegration.hpp"
#include "active_model.hpp"

using std::array;

void calculate_discretization(
    Model &model,
    double &sigma,
    Eigen::Matrix<double, Model::n_states, K> &X,
    Eigen::Matrix<double, Model::n_inputs, K> &U,
    array<Model::StateMatrix, (K - 1)> &A_bar,
    array<Model::ControlMatrix, (K - 1)> &B_bar,
    array<Model::ControlMatrix, (K - 1)> &C_bar,
    array<Model::StateVector, (K - 1)> &Sigma_bar,
    array<Model::StateVector, (K - 1)> &z_bar);