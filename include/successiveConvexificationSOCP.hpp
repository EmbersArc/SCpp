#pragma once

#include <array>

#include "optimizationProblem.hpp"
#include "active_model.hpp"

using std::array;

namespace sc
{
optimization_problem::SecondOrderConeProgram
build_successive_convexification_SOCP(
    Model &model,
    double &weight_trust_region_sigma,
    double &weight_trust_region_xu,
    double &weight_virtual_control,
    Eigen::Matrix<double, Model::state_dim_, K> &X,
    Eigen::Matrix<double, Model::input_dim_, K> &U,
    double &sigma,
    array<Model::state_matrix_t, (K - 1)> &A_bar,
    array<Model::control_matrix_t, (K - 1)> &B_bar,
    array<Model::control_matrix_t, (K - 1)> &C_bar,
    array<Model::state_vector_t, (K - 1)> &Sigma_bar,
    array<Model::state_vector_t, (K - 1)> &z_bar);
}