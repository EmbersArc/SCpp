#pragma once


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

    using StateVector   = Eigen::Matrix<double, n_states,        1>;
    using ControlVector = Eigen::Matrix<double, n_inputs,        1>;
    using StateMatrix   = Eigen::Matrix<double, n_states, n_states>;
    using ControlMatrix = Eigen::Matrix<double, n_states, n_inputs>;

    double total_time_guess() { return 3; }

    template<int K>
    void initialize(Eigen::Matrix<double, n_states, K> &X, Eigen::Matrix<double, n_inputs, K> &U) {

        //initial state
        double m_wet = 2;
        ControlVector r_I_init(4., 4., 0.);
        ControlVector v_I_init(0., -2., -2.);
        Vector4d q_B_I_init(1.0, 0.0, 0.0, 0.0);
        ControlVector w_B_init(0., 0., 0.);
        VectorXd x_init(14);

        //final state
        double m_dry = 1;
        ControlVector r_I_final(0., 0., 0.);
        ControlVector v_I_final(-1e-1, 0., 0.);
        Vector4d q_B_I_final(1.0, 0.0, 0.0, 0.0);
        ControlVector w_B_final(0., 0., 0.);
        VectorXd x_final(14);


        //gravity vector
        ControlVector g_I(-1, 0, 0);

        x_init << m_wet, r_I_init, v_I_init, q_B_I_init, w_B_init;
        x_final << m_dry, r_I_final, v_I_final, q_B_I_final, w_B_final;



        double alpha1, alpha2;
        for(int k=0; k<K; k++) {
            alpha1 = double(K-k)/K;
            alpha2 = double(k)/K;
            X(0, k) = alpha1 * x_init(0) + alpha2 * x_final(0);
            X.col(k).segment(1, 6) = alpha1 * x_init.segment(1, 6) + alpha2 * x_final.segment(1, 6);
            X.col(k).segment(7, 4) << 1., 0., 0., 0.;
            X.col(k).segment(11, 3) = alpha1 * x_init.segment(11, 3) + alpha2 * x_final.segment(11, 3);

            U.col(k) = X(0, k) * -g_I;
        }
    }

    StateVector                ode(const StateVector &x, const ControlVector &u);
    StateMatrix     state_jacobian(const StateVector &x, const ControlVector &u);
    ControlMatrix control_jacobian(const StateVector &x, const ControlVector &u);

    void add_application_constraints(EcosWrapper &solver, size_t K);

    StateVector get_random_state();
    ControlVector get_random_input();

};