#pragma once

#include "constants.hpp"
#include "EcosWrapper.hpp"
#include <Eigen/Dense>
using namespace Eigen;


class model_landing_6dof {

    const Vector3d g_I = Vector3d(-1, 0, 0);
    const Vector3d J_B = Vector3d(1e-2, 1e-2, 1e-2);
    const Vector3d r_T_B = Vector3d(-1e-2, 0, 0);
    const double alpha_m = 0.02;

public:

    static constexpr size_t n_states = 14;
    static constexpr size_t n_inputs = 3;
    static string get_name(){ return "model_landing_6dof"; }

    using StateVector   = Eigen::Matrix<double, n_states,        1>;
    using ControlVector = Eigen::Matrix<double, n_inputs,        1>;
    using StateMatrix   = Eigen::Matrix<double, n_states, n_states>;
    using ControlMatrix = Eigen::Matrix<double, n_states, n_inputs>;

    double total_time_guess() { return 3; }

    void initialize(Eigen::Matrix<double, n_states, K> &X, Eigen::Matrix<double, n_inputs, K> &U);

    StateVector                ode(const StateVector &x, const ControlVector &u);
    StateMatrix     state_jacobian(const StateVector &x, const ControlVector &u);
    ControlMatrix control_jacobian(const StateVector &x, const ControlVector &u);

    void add_application_constraints(optimization_problem::SecondOrderConeProgram &socp);

    StateVector get_random_state();
    ControlVector get_random_input();

};