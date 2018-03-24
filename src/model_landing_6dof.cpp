#include "model_landing_6dof.h"



void model_landing_6dof::initialize(int K, MatrixXd &X, MatrixXd &U) {

    //initial state
    double m_wet = 2;
    Vector3d r_I_init(4., 4., 0.);
    Vector3d v_I_init(0., -2., -2.);
    Vector4d q_B_I_init(1.0, 0.0, 0.0, 0.0);
    Vector3d w_B_init(0., 0., 0.);
    VectorXd x_init(14);

    //final state
    double m_dry = 1;
    Vector3d r_I_final(0., 0., 0.);
    Vector3d v_I_final(-1e-1, 0., 0.);
    Vector4d q_B_I_final(1.0, 0.0, 0.0, 0.0);
    Vector3d w_B_final(0., 0., 0.);
    VectorXd x_final(14);


    //gravity vector
    Vector3d g_I(-1, 0, 0);
    
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