#pragma once

#include "activeModelSC.hpp"

// void exactLinearDiscretization(Model &model,
//                                double &ts,
//                                const Model::state_vector_t &x_eq,
//                                const Model::input_vector_t &u_eq,
//                                Model::state_matrix_t &A,
//                                Model::control_matrix_t &B);
/**
 * @brief Multiple-shooting discretization with first order hold on input
 * 
 * @param model 
 * @param sigma 
 * @param X 
 * @param U 
 * @param A_bar 
 * @param B_bar 
 * @param C_bar 
 * @param S_bar 
 * @param z_bar 
 */
void calculateDiscretization(
    Model &model,
    double &sigma,
    const Eigen::MatrixXd &X,
    const Eigen::MatrixXd &U,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &S_bar,
    Model::state_vector_v_t &z_bar);
