#pragma once

#include "OptimizationProblem.hpp"
#include "active_model.hpp"
#include <array>

using std::array;

optimization_problem::SecondOrderConeProgram build_successive_convexification_SOCP (
    Model &model,
    double &weight_trust_region_sigma,
    double &weight_trust_region_xu,
    double &weight_virtual_control,
    Eigen::Matrix<double, Model::n_states, K> &X,
    Eigen::Matrix<double, Model::n_inputs, K> &U,
    double &sigma,
    array<Model::StateMatrix,   (K-1)> &A_bar,
    array<Model::ControlMatrix, (K-1)> &B_bar,
    array<Model::ControlMatrix, (K-1)> &C_bar,
    array<Model::StateVector,   (K-1)> &Sigma_bar,
    array<Model::StateVector,   (K-1)> &z_bar
);