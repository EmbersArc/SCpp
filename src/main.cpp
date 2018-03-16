#include <params.h>

using namespace std;
using namespace boost::numeric::odeint;

runge_kutta_dopri5<MatrixXd, double, MatrixXd, double, vector_space_algebra> stepper;


class ode_dPhidt{
    double sigma;
    Vector3d u_t, u_t1;
public:
    ode_dPhidt(const Vector3d &u_t, const Vector3d &u_t1, const double sigma) : u_t(u_t), u_t1(u_t1), sigma(sigma){}

    void operator()(const MatrixXd &V, MatrixXd &dVdt, const double t){
        Vector3d u = u_t + t / dt * (u_t1 - u_t);
        cout << dVdt.cols() << dVdt.rows() << endl;  //prints 00 instead of 1415
        dVdt.col(0) = sigma * f(V.col(0), u);
        dVdt.rightCols(14) = A(V.col(0), u, sigma) * V.rightCols(14);

    }
};

class ode_dVdt{
    double sigma;
    Vector3d u_t, u_t1;
public:
    ode_dVdt(const Vector3d &u_t, const Vector3d &u_t1, const double sigma) : u_t(u_t), u_t1(u_t1), sigma(sigma){}

    void operator()(const MatrixXd &V, MatrixXd &dVdt, const double t){

        Vector3d u = u_t + t / dt * (u_t1 - u_t);

        double alpha = t / dt;
        double beta = 1 - alpha;

        MatrixXd Phi_A_xi(14,15);
        Phi_A_xi.col(0) = V.col(0);
        Phi_A_xi.rightCols(14).setIdentity();

        //ode_dPhidt dPhidt(u_t, u_t1, sigma);
        //integrate_const(stepper, dPhidt, Phi_A_xi, t, dt, dt/10.);
        
        cout << dVdt.cols() << dVdt.rows() << endl;  //prints 00 instead of 1415

        dVdt.block<14, 1>(0, 0) = sigma * f(V.col(0), u);
        dVdt.block<14, 14>(0, 1) = A(V.col(0), u, sigma) * Phi_A_xi.block<14,14>(0,1);
        dVdt.block<14, 3>(0, 15) = Phi_A_xi.block<14,14>(0,1) * B(V.col(0), u, sigma) * alpha;
        dVdt.block<14, 3>(0, 18) = Phi_A_xi.block<14,14>(0,1) * B(V.col(0), u, sigma) * beta;
        dVdt.block<14, 1>(0, 21) = Phi_A_xi.block<14,14>(0,1) * f(V.col(0), u);
        dVdt.block<14, 1>(0, 22) = Phi_A_xi.block<14,14>(0,1) * (-A(V.col(0), u, sigma) * V.col(0) - B(V.col(0), u, sigma) * u);
    }
};

int main(){

    x_init << m_wet, r_I_init, v_I_init, q_B_I_init, w_B_init;
    x_final << m_dry, r_I_final, v_I_final, q_B_I_final, w_B_final;

    MatrixXd X(K, 14);
    MatrixXd U(K, 3);

//START INITIALIZATION
    cout << "Starting initialization." << endl;
    double alpha1, alpha2;
    for(int k=0; k<K; k++){
        alpha1 = double(K-k)/K;
        alpha2 = double(k)/K;
        X(k, 0) = alpha1 * x_init(0) + alpha2 * x_final(0);
        X.row(k).segment(1, 6) = alpha1 * x_init.segment(1, 6) + alpha2 * x_final.segment(1, 6);
        X.row(k).segment(7, 4) << 1, 0, 0, 0;
        X.row(k).segment(11, 3) = alpha1 * x_init.segment(11, 3) + alpha2 * x_final.segment(11, 3);

        U.row(k) = X(k, 0) * -g_I;
    }

    cout << "Initialization finished." << endl;


//START SUCCESSIVE CONVEXIFICATION
    double sigma = sigma_guess;

    for(int it = 1; it < iterations + 1; it++) {
        cout << "Iteration " << it << endl;
        vector<Matrix14d> A_bar(K);
        vector<Matrix14x3d> B_bar(K);
        vector<Matrix14x3d> C_bar(K);
        vector<Vector14d> Sigma_bar(K);
        vector<Vector14d> z_bar(K);

        cout << "Calculating new transition matrices." << endl;
        Matrix<double, 14, 23> V;

        for (int k = 0; k < K-1; k++) {
            V.setZero();
            V.col(0) = X.row(k);
            V.block<14,14>(0, 1).setIdentity();

            ode_dVdt dVdt(U.row(k), U.row(k+1), sigma);
            integrate_const(stepper, dVdt, V, 0., dt, dt/10.);

        }
    }

}