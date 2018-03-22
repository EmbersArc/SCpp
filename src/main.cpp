#include <params.h>


AMatrix A;
BMatrix B;
fVector f;

typedef Eigen::Matrix<double,14,23> state_type;
class ode_dVdt{
private:
    Vector3d u_t, u_t1;
    double sigma;
    Matrix14d Phi_A_xi;
    double alpha, beta;

public:
    void Update(const Vector3d &u_t, const Vector3d &u_t1, const double &sigma){
        this->u_t = u_t;
        this->u_t1 = u_t1;
        this->sigma = sigma;
    }

    void operator()(const state_type &V, state_type &dVdt, const double t){
        dVdt.setZero();
        Vector3d u = u_t + t / dt * (u_t1 - u_t);
        A.Update(V.col(0), u, sigma);
        B.Update(V.col(0), u, sigma);
        f.Update(V.col(0), u);

        alpha = t / dt;
        beta = 1. - alpha;

        Phi_A_xi = V.block<14, 14>(0, 1).inverse();

        dVdt.block<14, 1>(0, 0) = sigma * f();
        dVdt.block<14, 14>(0, 1) = A()* Phi_A_xi;
        dVdt.block<14, 3>(0, 15) = Phi_A_xi * B() * alpha;
        dVdt.block<14, 3>(0, 18) = Phi_A_xi * B() * beta;
        dVdt.block<14, 1>(0, 21) = Phi_A_xi * f();
        dVdt.block<14, 1>(0, 22) = Phi_A_xi * (-A() * V.col(0) - B() * u);

    }
};

int main(){

    x_init << m_wet, r_I_init, v_I_init, q_B_I_init, w_B_init;
    x_final << m_dry, r_I_final, v_I_final, q_B_I_final, w_B_final;

    MatrixXd X(14, K);
    MatrixXd U(3, K);


//START INITIALIZATION
    cout << "Starting initialization." << endl;
    double alpha1, alpha2;
    for(int k=0; k<K; k++){
        alpha1 = double(K-k)/K;
        alpha2 = double(k)/K;
        X(0, k) = alpha1 * x_init(0) + alpha2 * x_final(0);
        X.col(k).segment(1, 6) = alpha1 * x_init.segment(1, 6) + alpha2 * x_final.segment(1, 6);
        X.col(k).segment(7, 4) << 1., 0., 0., 0.;
        X.col(k).segment(11, 3) = alpha1 * x_init.segment(11, 3) + alpha2 * x_final.segment(11, 3);

        U.col(k) = X(0, k) * -g_I;
    }

    cout << "Initialization finished." << endl;

//START SUCCESSIVE CONVEXIFICATION
    double sigma;
    sigma = sigma_guess;

    state_type V;

    vector<Matrix14d> A_bar(K);
    vector<Matrix14x3d> B_bar(K);
    vector<Matrix14x3d> C_bar(K);
    vector<Vector14d> Sigma_bar(K);
    vector<Vector14d> z_bar(K);

    runge_kutta_dopri5<state_type, double, state_type, double, vector_space_algebra> stepper;
    ode_dVdt dVdt;

    for(int it = 1; it < iterations + 1; it++) {
        cout << "Iteration " << it << endl;
        cout << "Calculating new transition matrices." << endl;

        const clock_t begin_time = clock();

        for (int k = 0; k < K-1; k++) {
            V.setZero();
            V.col(0) = X.col(k);
            V.block<14,14>(0, 1).setIdentity();

            dVdt.Update(U.col(k), U.col(k+1), sigma);
            integrate_adaptive(make_controlled(1E-12 , 1E-12 , stepper), dVdt, V, 0., dt, dt/10.);

            A_bar[k] = V.block<14,14>(0, 1) * V.block<14,14>(0, 1);
            B_bar[k] = V.block<14,14>(0, 1) * V.block<14,3>(0, 15);
            C_bar[k] = V.block<14,14>(0, 1) * V.block<14,3>(0, 18);
            Sigma_bar[k] = V.block<14,14>(0, 1) * V.block<14,1>(0, 21);
            z_bar[k] = V.block<14,14>(0, 1) * V.block<14,1>(0, 22);
        }
        cout << "Transition matrices calculated in " << double( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds." << endl;

        // TODO: Solve problem.


    }

}