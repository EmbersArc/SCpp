#pragma once

#include "constants.hpp"
#include "EcosWrapper.hpp"
#include <Eigen/Dense>
using namespace Eigen;

class model_landing_6dof {

public:
    static constexpr size_t n_states = 14;
    static constexpr size_t n_inputs = 3;

    using StateVector   = Eigen::Matrix<double, n_states,        1>;
    using ControlVector = Eigen::Matrix<double, n_inputs,        1>;
    using StateMatrix   = Eigen::Matrix<double, n_states, n_states>;
    using ControlMatrix = Eigen::Matrix<double, n_states, n_inputs>;

    static string get_name(){ return "model_landing_6dof"; }

    double total_time_guess() { return 7; }

    void initialize(Eigen::Matrix<double, n_states, K> &X, Eigen::Matrix<double, n_inputs, K> &U);

    StateVector                ode(const StateVector &x, const ControlVector &u);
    StateMatrix     state_jacobian(const StateVector &x, const ControlVector &u);
    ControlMatrix control_jacobian(const StateVector &x, const ControlVector &u);

    void add_application_constraints(optimization_problem::SecondOrderConeProgram &socp,
                                     Eigen::Matrix<double, n_states, K> &X0,
                                     Eigen::Matrix<double, n_inputs, K> &U0);

    StateVector get_random_state();
    ControlVector get_random_input();

private:
    const Vector3d g_I = Vector3d(-1, 0, 0);
    const Vector3d J_B = Vector3d(1e-2, 1e-2, 1e-2);
    const Vector3d r_T_B = Vector3d(-1e-2, 0, 0);
    const double alpha_m = 0.1;
    const double T_min = 0.3;
    const double T_max = 5.;

    //initial state
    const double m_wet = 2.;
    Vector3d r_I_init = Vector3d(10, 5, 2);
    Vector3d v_I_init = Vector3d(-5, -4, -2);
    Vector4d q_B_I_init = Vector4d(1.0, 0.0, 0.0, 0.0);
    Vector3d w_B_init = Vector3d(0., 0., 0.);

    //final state
    const double m_dry = 1.;
    Vector3d r_I_final = Vector3d(0., 0., 0.);
    Vector3d v_I_final = Vector3d(-1e-1, 0., 0.);
    Vector4d q_B_I_final = Vector4d(1.0, 0.0, 0.0, 0.0);
    Vector3d w_B_final = Vector3d(0., 0., 0.);

    const double tan_delta_max = tan(20. / 180. * PI);
    const double cos_theta_max = cos(90. / 180. * PI);
    const double tan_gamma_gs = tan(20. / 180. * PI);
    const double w_B_max = 60. / 180. * PI;

};