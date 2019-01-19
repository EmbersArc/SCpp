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
    Eigen::MatrixXd &X,
    Eigen::MatrixXd &U,
    double &sigma,
    vector<Model::state_matrix_t> &A_bar,
    vector<Model::control_matrix_t> &B_bar,
    vector<Model::control_matrix_t> &C_bar,
    vector<Model::state_vector_t> &Sigma_bar,
    vector<Model::state_vector_t> &z_bar);
}