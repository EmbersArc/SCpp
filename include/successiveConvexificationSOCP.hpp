#pragma once

#include "optimizationProblem.hpp"
#include "active_model.hpp"

using std::array;

namespace sc
{
op::SecondOrderConeProgram build_sc_SOCP(
    Model &model,
    double &weight_trust_region_sigma,
    Model::state_vector_t &weight_trust_region_x,
    Model::input_vector_t &weight_trust_region_u,
    Model::state_vector_t &weight_virtual_control,
    Eigen::MatrixXd &X,
    Eigen::MatrixXd &U,
    double &sigma,
    vector<Model::state_matrix_t> &A_bar,
    vector<Model::control_matrix_t> &B_bar,
    vector<Model::control_matrix_t> &C_bar,
    vector<Model::state_vector_t> &Sigma_bar,
    vector<Model::state_vector_t> &z_bar);
}