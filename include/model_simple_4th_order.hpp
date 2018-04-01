#pragma once

#include "constants.hpp"
#include "EcosWrapper.hpp"
#include <Eigen/Dense>


class model_simple_4th_order {
public:

    static constexpr size_t n_states = 4;
    static constexpr size_t n_inputs = 1;
    static string get_name(){ return "model_simple_4th_order"; }

    using StateVector   = Eigen::Matrix<double, n_states,        1>;
    using ControlVector = Eigen::Matrix<double, n_inputs,        1>;
    using StateMatrix   = Eigen::Matrix<double, n_states, n_states>;
    using ControlMatrix = Eigen::Matrix<double, n_states, n_inputs>;

    double total_time_guess() { return 3; }

    void initialize(Eigen::Matrix<double, n_states, K> &X, Eigen::Matrix<double, n_inputs, K> &U) {
        X.setZero();
        U.setZero();
    }

    StateVector                ode(const StateVector &x, const ControlVector &u);
    StateMatrix     state_jacobian(const StateVector &x, const ControlVector &u);
    ControlMatrix control_jacobian(const StateVector &x, const ControlVector &u);
    
    void add_application_constraints(optimization_problem::SecondOrderConeProgram &socp);
    
    StateVector get_random_state();
    ControlVector get_random_input();

};