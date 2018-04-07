#include "model_landing_6dof.h"



void model_landing_6dof::initialize(Eigen::Matrix<double, n_states, K> &X, Eigen::Matrix<double, n_inputs, K> &U) {

    double alpha1, alpha2;

    x_init << m_wet, r_I_init, v_I_init, q_B_I_init, w_B_init;
    x_final << m_dry, r_I_final, v_I_final, q_B_I_final, w_B_final;

    for(int k=0; k<K; k++) {
        alpha1 = double(K - k) / K;
        alpha2 = double(k) / K;
        X(0, k) = alpha1 * x_init(0) + alpha2 * x_final(0);
        X.col(k).segment(1, 6) = alpha1 * x_init.segment(1, 6) + alpha2 * x_final.segment(1, 6);
        X.col(k).segment(7, 4) << 1., 0., 0., 0.;
        X.col(k).segment(11, 3) = alpha1 * x_init.segment(11, 3) + alpha2 * x_final.segment(11, 3);

        U.col(k) = X(0, k) * -g_I;
    }
}


model_landing_6dof::StateMatrix model_landing_6dof::state_jacobian(const StateVector &x, const ControlVector &u) {
    StateMatrix A;

    double t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24;
    double t25, t26, t27, t28, t29, t30, t31, t32;

    A.setZero();
    t2 = 1.0 / (x[0] * x[0]);
    t3 = 1.0 / x[0];
    t4 = x[10] * x[10];
    t5 = t4 * 2.0;
    t6 = x[7] * x[10] * 2.0;
    t7 = x[8] * x[9] * 2.0;
    t8 = t3 * u[2] * x[8] * 2.0;
    t9 = t3 * u[2] * x[10] * 2.0;
    t10 = x[8] * x[8];
    t11 = t10 * 2.0;
    t12 = x[9] * x[9];
    t13 = t12 * 2.0;
    t14 = x[7] * x[9] * 2.0;
    t15 = x[7] * x[8] * 2.0;
    t16 = x[9] * x[10] * 2.0;
    t17 = t3 * u[0] * x[9] * 2.0;
    t18 = t3 * u[1] * x[8] * 2.0;
    t19 = t3 * u[1] * x[7] * 2.0;
    t20 = t3 * u[1] * x[10] * 2.0;
    t21 = t3 * u[0] * x[8] * 2.0;
    t22 = t3 * u[1] * x[9] * 2.0;
    t23 = x[13] * (1.0 / 2.0);
    t24 = x[11] * (1.0 / 2.0);
    t25 = x[7] * (1.0 / 2.0);
    t26 = x[12] * (1.0 / 2.0);
    t27 = x[9] * (1.0 / 2.0);
    t28 = 1.0 / J_B[0];
    t29 = 1.0 / J_B[1];
    t30 = J_B[1] * x[12];
    t31 = 1.0 / J_B[2];
    t32 = J_B[0] * x[11];
    A(1, 4) = 1.0;
    A(2, 5) = 1.0;
    A(3, 6) = 1.0;
    A(4, 0) = t2 * u[2] * (t14 - x[8] * x[10] * 2.0) - t2 * u[1] * (t6 + t7) + t2 * u[0] * (t5 + t13 - 1.0);
    A(4, 7) = t20 - t3 * u[2] * x[9] * 2.0;
    A(4, 8) = t9 + t22;
    A(4, 9) = t18 - t3 * u[0] * x[9] * 4.0 - t3 * u[2] * x[7] * 2.0;
    A(4, 10) = t8 + t19 - t3 * u[0] * x[10] * 4.0;
    A(5, 0) = t2 * u[1] * (t5 + t11 - 1.0) - t2 * u[2] * (t15 + t16) + t2 * u[0] * (t6 - t7);
    A(5, 7) = t8 - t3 * u[0] * x[10] * 2.0;
    A(5, 8) = t17 - t3 * u[1] * x[8] * 4.0 + t3 * u[2] * x[7] * 2.0;
    A(5, 9) = t9 + t21;
    A(5, 10) = t3 * u[0] * x[7] * -2.0 - t3 * u[1] * x[10] * 4.0 + t3 * u[2] * x[9] * 2.0;
    A(6, 0) = -t2 * u[0] * (t14 + x[8] * x[10] * 2.0) + t2 * u[2] * (t11 + t13 - 1.0) + t2 * u[1] * (t15 - t16);
    A(6, 7) = t17 - t18;
    A(6, 8) = -t19 + t3 * u[0] * x[10] * 2.0 - t3 * u[2] * x[8] * 4.0;
    A(6, 9) = t20 + t3 * u[0] * x[7] * 2.0 - t3 * u[2] * x[9] * 4.0;
    A(6, 10) = t21 + t22;
    A(7, 8) = x[11] * (-1.0 / 2.0);
    A(7, 9) = x[12] * (-1.0 / 2.0);
    A(7, 10) = x[13] * (-1.0 / 2.0);
    A(7, 11) = x[8] * (-1.0 / 2.0);
    A(7, 12) = x[9] * (-1.0 / 2.0);
    A(7, 13) = x[10] * (-1.0 / 2.0);
    A(8, 7) = t24;
    A(8, 9) = t23;
    A(8, 10) = x[12] * (-1.0 / 2.0);
    A(8, 11) = t25;
    A(8, 12) = x[10] * (-1.0 / 2.0);
    A(8, 13) = t27;
    A(9, 7) = t26;
    A(9, 8) = -t23;
    A(9, 10) = t24;
    A(9, 11) = x[10] * (1.0 / 2.0);
    A(9, 12) = t25;
    A(9, 13) = x[8] * (-1.0 / 2.0);
    A(10, 7) = t23;
    A(10, 8) = t26;
    A(10, 9) = -t24;
    A(10, 11) = -t27;
    A(10, 12) = x[8] * (1.0 / 2.0);
    A(10, 13) = t25;
    A(11, 12) = t28 * (J_B[1] * x[13] - J_B[2] * x[13]);
    A(11, 13) = t28 * (t30 - J_B[2] * x[12]);
    A(12, 11) = -t29 * (J_B[0] * x[13] - J_B[2] * x[13]);
    A(12, 13) = -t29 * (t32 - J_B[2] * x[11]);
    A(13, 11) = -t31 * (t30 - J_B[0] * x[12]);
    A(13, 12) = t31 * (t32 - J_B[1] * x[11]);

    return A;
}



model_landing_6dof::ControlMatrix model_landing_6dof::control_jacobian(const StateVector &x, const ControlVector &u) {
    ControlMatrix B;
    B.setZero();
    B.row(0) << -(alpha_m * u[0]) / u.norm(), -(alpha_m * u[1]) / u.norm(), -(alpha_m * u[2]) / u.norm();
    B.row(4) << -(2 * pow(x[9], 2) + 2 * pow(x[10], 2) - 1) / x[0],
            (2 * x[7] * x[10] + 2 * x[8] * x[9]) / x[0],
            -(2 * x[7] * x[9] - 2 * x[8] * x[10]) / x[0];
    B.row(5) << -(2 * x[7] * x[10] - 2 * x[8] * x[9]) / x[0],
            -(2 * pow(x[8], 2) + 2 * pow(x[10], 2) - 1) / x[0],
            (2 * x[7] * x[8] + 2 * x[9] * x[10]) / x[0];
    B.row(6) << (2 * x[7] * x[9] + 2 * x[8] * x[10]) / x[0],
            -(2 * x[7] * x[8] - 2 * x[9] * x[10]) / x[0],
            -(2 * pow(x[8], 2) + 2 * pow(x[9], 2) - 1) / x[0];
    B.row(11) << 0., -r_T_B[2] / J_B[0], r_T_B[1] / J_B[0];
    B.row(12) << r_T_B[2] / J_B[1], 0., -r_T_B[0] / J_B[1];
    B.row(13) << -r_T_B[1] / J_B[2], r_T_B[0] / J_B[2], 0.;
    return B; 
}



model_landing_6dof::StateVector model_landing_6dof::ode(const StateVector &x, const ControlVector &u) {
    StateVector f;
    f <<
            -alpha_m * u.norm(),
            x[4],
            x[5],
            x[6],
            g_I[0] - (u[0] * (2 * pow(x[9], 2) + 2 * pow(x[10], 2) - 1)) / x[0] + (u[1] * (2 * x[7] * x[10] + 2 * x[8] * x[9])) / x[0] - (u[2] * (2 * x[7] * x[9] - 2 * x[8] * x[10])) / x[0],
            g_I[1] - (u[1] * (2 * pow(x[8], 2) + 2 * pow(x[10], 2) - 1)) / x[0] - (u[0] * (2 * x[7] * x[10] - 2 * x[8] * x[9])) / x[0] + (u[2] * (2 * x[7] * x[8] + 2 * x[9] * x[10])) / x[0],
            g_I[2] - (u[2] * (2 * pow(x[8], 2) + 2 * pow(x[9], 2) - 1)) / x[0] + (u[0] * (2 * x[7] * x[9] + 2 * x[8] * x[10])) / x[0] - (u[1] * (2 * x[7] * x[8] - 2 * x[9] * x[10])) / x[0],
            - (x[8] * x[11]) / 2 - (x[9] * x[12]) / 2 - (x[10] * x[13]) / 2,
            (x[7] * x[11]) / 2 + (x[9] * x[13]) / 2 - (x[10] * x[12]) / 2,
            (x[7] * x[12]) / 2 - (x[8] * x[13]) / 2 + (x[10] * x[11]) / 2,
            (x[7] * x[13]) / 2 + (x[8] * x[12]) / 2 - (x[9] * x[11]) / 2,
            (r_T_B[1] * u[2] - r_T_B[2] * u[1] + J_B[1] * x[12] * x[13] - J_B[2] * x[12] * x[13]) / J_B[0],
            - (r_T_B[0] * u[2] - r_T_B[2] * u[0] + J_B[0] * x[11] * x[13] - J_B[2] * x[11] * x[13]) / J_B[1],
            (r_T_B[0] * u[1] - r_T_B[1] * u[0] + J_B[0] * x[11] * x[12] - J_B[1] * x[11] * x[12]) / J_B[2];
    return f;
}

void model_landing_6dof::add_application_constraints(
        optimization_problem::SecondOrderConeProgram &socp,
        Eigen::Matrix<double, n_states, K> &X0,
        Eigen::Matrix<double, n_inputs, K> &U0
) {
    auto var = [&](const string &name, const vector<size_t> &indices){ return socp.get_variable(name,indices); };
    auto param = [](double &param_value){ return optimization_problem::Parameter(&param_value); };
    auto param_fn = [](std::function<double()> callback){ return optimization_problem::Parameter(callback); };

    x_init << m_wet, r_I_init, v_I_init, q_B_I_init, w_B_init;
    x_final << m_dry, r_I_final, v_I_final, q_B_I_final, w_B_final;

    // initial state
    socp.add_constraint( (-1.0) * var("X", {0, 0}) + (x_init(0)) == 0.0 );
    socp.add_constraint( (-1.0) * var("X", {1, 0}) + (x_init(1)) == 0.0 );
    socp.add_constraint( (-1.0) * var("X", {2, 0}) + (x_init(2)) == 0.0 );
    socp.add_constraint( (-1.0) * var("X", {3, 0}) + (x_init(3)) == 0.0 );
    socp.add_constraint( (-1.0) * var("X", {4, 0}) + (x_init(4)) == 0.0 );
    socp.add_constraint( (-1.0) * var("X", {5, 0}) + (x_init(5)) == 0.0 );
    socp.add_constraint( (-1.0) * var("X", {6, 0}) + (x_init(6)) == 0.0 );
//    socp.add_constraint( (-1.0) * var("X", {7, 0}) + (x_init(7)) == 0.0 );
//    socp.add_constraint( (-1.0) * var("X", {8, 0}) + (x_init(8)) == 0.0 );
//    socp.add_constraint( (-1.0) * var("X", {9, 0}) + (x_init(9)) == 0.0 );
//    socp.add_constraint( (-1.0) * var("X", {10, 0}) + (x_init(10)) == 0.0 );
    socp.add_constraint( (-1.0) * var("X", {11, 0}) + (x_init(11)) == 0.0 );
    socp.add_constraint( (-1.0) * var("X", {12, 0}) + (x_init(12)) == 0.0 );
    socp.add_constraint( (-1.0) * var("X", {13, 0}) + (x_init(13)) == 0.0 );


    // Final State
    for(size_t i = 0; i<n_states; i++){
        socp.add_constraint( (-1.0) * var("X", {i, K-1}) + (x_final(i)) == 0.0 );
    }
    socp.add_constraint( (1.0) * var("U", {1, K-1}) == (0.0) );
    socp.add_constraint( (1.0) * var("U", {2, K-1}) == (0.0) );


    // State Constraints:
    for(size_t k = 0; k<K; k++){

        // Mass
        //     x(0) >= m_dry
        //     for all k
        socp.add_constraint( (1.0) * var("X", {0, k}) + (-m_dry) >= (0.0) );

        // Max Tilt Angle
        //
        // x(9) ^ 2 + x(10) ^ 2 <= c
        // with c := -(cos_theta_max - 1) / 2
        //
        // equivalent to
        // norm2(
        //          (1 + c) / 2
        //              x(9)
        //              x(10)
        // )
        // <= (1 - c) / 2

        /*double c = -(cos_theta_max - 1.) / 2.;
        double c1 = (1. + c) / 2.;
        double c2 = (1. - c) / 2.;

        socp.add_constraint( optimization_problem::norm2({
                param(c1),
                (1.0) * var("X", {9, k}),
                (1.0) * var("X", {10, k})
        }) <= (c2) );*/

        // Glide Slope
        /*socp.add_constraint(
                optimization_problem::norm2({ (1.0) * var("X", {2, k}),
                                              (1.0) * var("X", {3, k}) })
                <= (1.0 / tan_gamma_gs) * var("X", {1, k})
        );

        // Max Rotation Velocity
        socp.add_constraint(
                optimization_problem::norm2({ (1.0) * var("X", {11, k}),
                                              (1.0) * var("X", {12, k}),
                                              (1.0) * var("X", {13, k}) })
                <= (w_B_max)
        );*/
    }


    // Control Constraints
    for(size_t k = 0; k<K; k++) {

        // Minimum Thrust
        /*optimization_problem::AffineExpression lhs;

        auto U0_norm = param_fn([&U0,k](){ return U0.norm(); });

        for (size_t i = 0; i < n_inputs; i++) {
            lhs = lhs + param(U0(i,k)) * var("U", {i, k});
        }

        socp.add_constraint(lhs + (-T_min) >= (0.0));*/


        // Maximum Thrust
        socp.add_constraint(
                optimization_problem::norm2({ (1.0) * var("U", {0, k}),
                                              (1.0) * var("U", {1, k}),
                                              (1.0) * var("U", {2, k}) })
                <= (T_max)
        );

        // Maximum Gimbal Angle
        socp.add_constraint(
                optimization_problem::norm2({ (1.0) * var("U", {0, k}),
                                              (1.0) * var("U", {1, k}),
                                              (1.0) * var("U", {2, k}) })
                <= (1.0 / cos_delta_max) * var("U", {0, k})
        );
    }
}

model_landing_6dof::StateVector model_landing_6dof::get_random_state() {
    StateVector X;
    X.setRandom();

    X(0) = abs(X(0)) + 1.;
    X(1) = abs(X(1)) + 1.;
    X.segment(7, 4).normalize();

    return X;
}

model_landing_6dof::ControlVector model_landing_6dof::get_random_input() {
    ControlVector U;
    U.setRandom();
    U.normalize();
    U(0) = abs(U(0));

    return U;
}

