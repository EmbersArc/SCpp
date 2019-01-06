#include "model_landing_6dof.h"



void model_landing_6dof::initialize(Eigen::Matrix<double, n_states, K> &X, Eigen::Matrix<double, n_inputs, K> &U) {

//    Nondimensionalize();

    StateVector x_init; x_init << m_wet, r_I_init, v_I_init, q_B_I_init, w_B_init;
    StateVector x_final; x_final << m_dry, r_I_final, v_I_final, q_B_I_final, w_B_final;

    for(int k=0; k<K; k++) {
        double alpha1 = double(K - k) / K;
        double alpha2 = double(k) / K;
        X(0, k) = alpha1 * x_init(0) + alpha2 * x_final(0);
        X.col(k).segment(1, 6) = alpha1 * x_init.segment(1, 6) + alpha2 * x_final.segment(1, 6);
        X.col(k).segment(7, 4) << 1., 0., 0., 0.;
        X.col(k).segment(11, 3) = alpha1 * x_init.segment(11, 3) + alpha2 * x_final.segment(11, 3);

        U.col(k) = X(0, k) * -g_I;
    }
}


model_landing_6dof::StateMatrix model_landing_6dof::state_jacobian(const StateVector &x, const ControlVector &u) {
    const double m = x(0, 0);
    const double q0 = x(7, 0);
    const double q1 = x(8, 0);
    const double q2 = x(9, 0);
    const double q3 = x(10, 0);
    const double wx = x(11, 0);
    const double wy = x(12, 0);
    const double wz = x(13, 0);
    const double ux = u(0, 0);
    const double uy = u(1, 0);
    const double uz = u(2, 0);
    const double J_B0 = J_B(0, 0);
    const double J_B1 = J_B(1, 0);
    const double J_B2 = J_B(2, 0);
    const double x0 = pow(m, -2);
    const double x1 = uy*x0;
    const double x2 = 2*q0;
    const double x3 = q3*x2;
    const double x4 = 2*q1;
    const double x5 = q2*x4;
    const double x6 = uz*x0;
    const double x7 = q2*x2;
    const double x8 = q3*x4;
    const double x9 = ux*x0;
    const double x10 = -2*pow(q2, 2);
    const double x11 = -2*pow(q3, 2);
    const double x12 = 1.0/m;
    const double x13 = 2*q2*x12;
    const double x14 = uz*x13;
    const double x15 = q3*uy*x12;
    const double x16 = 2*x15;
    const double x17 = uy*x13;
    const double x18 = uz*x12;
    const double x19 = 2*q3*x18;
    const double x20 = 2*q0*x12;
    const double x21 = uz*x20;
    const double x22 = 2*q1*x12;
    const double x23 = uy*x22;
    const double x24 = 4*q2;
    const double x25 = ux*x12;
    const double x26 = uy*x20;
    const double x27 = uz*x22;
    const double x28 = q3*ux*x12;
    const double x29 = q1*x2;
    const double x30 = 2*q2;
    const double x31 = q3*x30;
    const double x32 = -2*pow(q1, 2) + 1;
    const double x33 = 2*x28;
    const double x34 = 4*q1*x12;
    const double x35 = x25*x30;
    const double x36 = ux*x22;
    const double x37 = ux*x20;
    const double x38 = 0.5*wx;
    const double x39 = -x38;
    const double x40 = 0.5*wy;
    const double x41 = -x40;
    const double x42 = 0.5*wz;
    const double x43 = -x42;
    const double x44 = 0.5*q1;
    const double x45 = -x44;
    const double x46 = 0.5*q2;
    const double x47 = -x46;
    const double x48 = 0.5*q3;
    const double x49 = -x48;
    const double x50 = 0.5*q0;
    const double x51 = 1.0/J_B0;
    const double x52 = J_B2*wz;
    const double x53 = J_B1*wy;
    const double x54 = 1.0/J_B1;
    const double x55 = J_B0*wx;
    const double x56 = 1.0/J_B2;

    StateMatrix A;
    A.setZero();
    A(1, 4) = 1;
    A(2, 5) = 1;
    A(3, 6) = 1;
    A(4, 0) = -x1*(-x3 + x5) - x6*(x7 + x8) - x9*(x10 + x11 + 1);
    A(4, 7) = x14 - x16;
    A(4, 8) = x17 + x19;
    A(4, 9) = x21 + x23 - x24*x25;
    A(4, 10) = -x26 + x27 - 4*x28;
    A(5, 0) = -x1*(x11 + x32) - x6*(-x29 + x31) - x9*(x3 + x5);
    A(5, 7) = -x27 + x33;
    A(5, 8) = -uy*x34 - x21 + x35;
    A(5, 9) = x19 + x36;
    A(5, 10) = x14 - 4*x15 + x37;
    A(6, 0) = -x1*(x29 + x31) - x6*(x10 + x32) - x9*(-x7 + x8);
    A(6, 7) = x23 - x35;
    A(6, 8) = -uz*x34 + x26 + x33;
    A(6, 9) = x16 - x18*x24 - x37;
    A(6, 10) = x17 + x36;
    A(7, 8) = x39;
    A(7, 9) = x41;
    A(7, 10) = x43;
    A(7, 11) = x45;
    A(7, 12) = x47;
    A(7, 13) = x49;
    A(8, 7) = x38;
    A(8, 9) = x42;
    A(8, 10) = x41;
    A(8, 11) = x50;
    A(8, 12) = x49;
    A(8, 13) = x46;
    A(9, 7) = x40;
    A(9, 8) = x43;
    A(9, 10) = x38;
    A(9, 11) = x48;
    A(9, 12) = x50;
    A(9, 13) = x45;
    A(10, 7) = x42;
    A(10, 8) = x40;
    A(10, 9) = x39;
    A(10, 11) = x47;
    A(10, 12) = x44;
    A(10, 13) = x50;
    A(11, 12) = x51*(J_B1*wz - x52);
    A(11, 13) = x51*(-J_B2*wy + x53);
    A(12, 11) = x54*(-J_B0*wz + x52);
    A(12, 13) = x54*(J_B2*wx - x55);
    A(13, 11) = x56*(J_B0*wy - x53);
    A(13, 12) = x56*(-J_B1*wx + x55);
    return A;
}



model_landing_6dof::ControlMatrix model_landing_6dof::control_jacobian(const StateVector &x, const ControlVector &u) {
    const double m = x(0, 0);
    const double q0 = x(7, 0);
    const double q1 = x(8, 0);
    const double q2 = x(9, 0);
    const double q3 = x(10, 0);
    const double ux = u(0, 0);
    const double uy = u(1, 0);
    const double uz = u(2, 0);
    const double r_T_B0 = r_T_B(0, 0);
    const double r_T_B1 = r_T_B(1, 0);
    const double r_T_B2 = r_T_B(2, 0);
    const double J_B0 = J_B(0, 0);
    const double J_B1 = J_B(1, 0);
    const double J_B2 = J_B(2, 0);
    const double x0 = alpha_m/sqrt(pow(ux, 2) + pow(uy, 2) + pow(uz, 2));
    const double x1 = 1.0/m;
    const double x2 = -2*pow(q2, 2);
    const double x3 = -2*pow(q3, 2);
    const double x4 = 2*q0;
    const double x5 = q3*x4;
    const double x6 = 2*q1;
    const double x7 = q2*x6;
    const double x8 = q2*x4;
    const double x9 = q3*x6;
    const double x10 = -2*pow(q1, 2) + 1;
    const double x11 = q1*x4;
    const double x12 = 2*q2*q3;
    const double x13 = 1.0/J_B0;
    const double x14 = 1.0/J_B1;
    const double x15 = 1.0/J_B2;

    ControlMatrix B;
    B.setZero();
    B(0, 0) = -ux*x0;
    B(0, 1) = -uy*x0;
    B(0, 2) = -uz*x0;
    B(4, 0) = x1*(x2 + x3 + 1);
    B(4, 1) = x1*(-x5 + x7);
    B(4, 2) = x1*(x8 + x9);
    B(5, 0) = x1*(x5 + x7);
    B(5, 1) = x1*(x10 + x3);
    B(5, 2) = x1*(-x11 + x12);
    B(6, 0) = x1*(-x8 + x9);
    B(6, 1) = x1*(x11 + x12);
    B(6, 2) = x1*(x10 + x2);
    B(11, 1) = -r_T_B2*x13;
    B(11, 2) = r_T_B1*x13;
    B(12, 0) = r_T_B2*x14;
    B(12, 2) = -r_T_B0*x14;
    B(13, 0) = -r_T_B1*x15;
    B(13, 1) = r_T_B0*x15;
    return B;
}



model_landing_6dof::StateVector model_landing_6dof::ode(const StateVector &x, const ControlVector &u) {
    const double m = x(0, 0);
    const double vx = x(4, 0);
    const double vy = x(5, 0);
    const double vz = x(6, 0);
    const double q0 = x(7, 0);
    const double q1 = x(8, 0);
    const double q2 = x(9, 0);
    const double q3 = x(10, 0);
    const double wx = x(11, 0);
    const double wy = x(12, 0);
    const double wz = x(13, 0);
    const double ux = u(0, 0);
    const double uy = u(1, 0);
    const double uz = u(2, 0);
    const double g_I0 = g_I(0, 0);
    const double g_I1 = g_I(1, 0);
    const double g_I2 = g_I(2, 0);
    const double r_T_B0 = r_T_B(0, 0);
    const double r_T_B1 = r_T_B(1, 0);
    const double r_T_B2 = r_T_B(2, 0);
    const double J_B0 = J_B(0, 0);
    const double J_B1 = J_B(1, 0);
    const double J_B2 = J_B(2, 0);
    const double x0 = 1.0/m;
    const double x1 = uy*x0;
    const double x2 = 2*q0;
    const double x3 = q3*x2;
    const double x4 = 2*q1;
    const double x5 = q2*x4;
    const double x6 = uz*x0;
    const double x7 = q2*x2;
    const double x8 = q3*x4;
    const double x9 = ux*x0;
    const double x10 = -2*pow(q2, 2);
    const double x11 = -2*pow(q3, 2);
    const double x12 = q1*x2;
    const double x13 = 2*q2*q3;
    const double x14 = -2*pow(q1, 2) + 1;
    const double x15 = 0.5*q1;
    const double x16 = 0.5*q2;
    const double x17 = 0.5*q3;
    const double x18 = 0.5*q0;
    const double x19 = J_B1*wy;
    const double x20 = J_B2*wz;
    const double x21 = J_B0*wx;

    StateVector f;
    f.setZero();
    f(0, 0) = -alpha_m*sqrt(pow(ux, 2) + pow(uy, 2) + pow(uz, 2));
    f(1, 0) = vx;
    f(2, 0) = vy;
    f(3, 0) = vz;
    f(4, 0) = g_I0 + x1*(-x3 + x5) + x6*(x7 + x8) + x9*(x10 + x11 + 1);
    f(5, 0) = g_I1 + x1*(x11 + x14) + x6*(-x12 + x13) + x9*(x3 + x5);
    f(6, 0) = g_I2 + x1*(x12 + x13) + x6*(x10 + x14) + x9*(-x7 + x8);
    f(7, 0) = -wx*x15 - wy*x16 - wz*x17;
    f(8, 0) = wx*x18 - wy*x17 + wz*x16;
    f(9, 0) = wx*x17 + wy*x18 - wz*x15;
    f(10, 0) = -wx*x16 + wy*x15 + wz*x18;
    f(11, 0) = (r_T_B1*uz - r_T_B2*uy - wy*x20 + wz*x19)/J_B0;
    f(12, 0) = (-r_T_B0*uz + r_T_B2*ux + wx*x20 - wz*x21)/J_B1;
    f(13, 0) = (r_T_B0*uy - r_T_B1*ux - wx*x19 + wy*x21)/J_B2;
    return f;
}


void model_landing_6dof::add_application_constraints(
        optimization_problem::SecondOrderConeProgram &socp,
        Eigen::Matrix<double, n_states, K> &X0,
        Eigen::Matrix<double, n_inputs, K> &U0
) {
    auto var = [&](const string &name, const vector<size_t> &indices){ return socp.get_variable(name,indices); };
//    auto param = [](double &param_value){ return optimization_problem::Parameter(&param_value); };
//    auto param_fn = [](std::function<double()> callback){ return optimization_problem::Parameter(callback); };

    StateVector x_init; x_init << m_wet, r_I_init, v_I_init, q_B_I_init, w_B_init;
    StateVector x_final; x_final << m_dry, r_I_final, v_I_final, q_B_I_final, w_B_final;

    // Initial state
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


    // Final State (mass is free)
    for(size_t i = 1; i<n_states; i++){
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
        // norm2([x(9), x(10)]) <= c
        // with c := sqrt((1 - cos_theta_max) / 2)
        const double c = sqrt((1.0 - cos_theta_max) / 2.);
        socp.add_constraint( optimization_problem::norm2({ 
            (1.0) * var("X", {9, k}),  
            (1.0) * var("X", {10, k}) 
        }) <= (c) );

        // Glide Slope
        socp.add_constraint(
            optimization_problem::norm2({ 
                (1.0) * var("X", {2, k}),
                (1.0) * var("X", {3, k}) 
            })
            <= (1.0 / tan_gamma_gs) * var("X", {1, k})
        );

        // Max Rotation Velocity
        socp.add_constraint(
            optimization_problem::norm2({ 
                (1.0) * var("X", {11, k}),
                (1.0) * var("X", {12, k}),
                (1.0) * var("X", {13, k}) 
            })
            <= (w_B_max)
        );
    }


    // Control Constraints
    for(size_t k = 0; k<K; k++) {

        // Linearized Minimum Thrust
        /*optimization_problem::AffineExpression lhs;
        for (size_t i = 0; i < n_inputs; i++) {
            lhs = lhs + param_fn([&U0,i,k](){ return (U0(i,k) / sqrt(U0(0,k)*U0(0,k) + U0(1,k)*U0(1,k) + U0(2,k)*U0(2,k))  ); }) * var("U", {i, k});
        }
        socp.add_constraint(lhs + (-T_min) >= (0.0));*/


        // Simplified Minimum Thrust
        socp.add_constraint( (1.0) * var("U", {0, k}) + (-T_min) >= (0.0));


        // Maximum Thrust
        socp.add_constraint(
            optimization_problem::norm2({ 
                (1.0) * var("U", {0, k}),
                (1.0) * var("U", {1, k}),
                (1.0) * var("U", {2, k}) 
            })
            <= (T_max)
        );

        // Maximum Gimbal Angle
        socp.add_constraint(
            optimization_problem::norm2({
                (1.0) * var("U", {1, k}),
                (1.0) * var("U", {2, k}) 
            })
            <= (tan_delta_max) * var("U", {0, k})
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

void model_landing_6dof::Nondimensionalize() {

    double r_scale = r_I_init[0];
    double m_scale = m_wet;

    alpha_m *= r_scale;
    r_T_B /= r_scale;
    g_I /= r_scale;
    J_B /= (m_scale * r_scale * r_scale);

    m_wet /= m_scale;
    r_I_init /= r_scale;
    v_I_init /= r_scale;

    m_dry /= m_scale;
    r_I_final /= r_scale;
    v_I_final /= r_scale;

    T_max /= m_scale * r_scale;
    T_min /= m_scale * r_scale;

}
