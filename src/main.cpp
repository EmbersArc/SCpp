#include <params.h>

using namespace std;
using namespace boost::numeric::odeint;

typedef Matrix<double,14,23> state_type1;
typedef Matrix<double,14,15> state_type2;

AMatrix A;
BMatrix B;
fVector f;

class ode_dPhidt{
    double sigma;
    Vector3d u_t, u_t1;
public:
    ode_dPhidt(Vector3d u_t, Vector3d u_t1, double sigma) : u_t(u_t), u_t1(u_t1), sigma(sigma){}

    void operator()(const state_type2 &V, state_type2 &dVdt, const double t){
        Vector3d u = u_t + t / dt * (u_t1 - u_t);
        A.Update(V.col(0), u, sigma);
        f.Update(V.col(0), u);

        dVdt.col(0) = sigma * f();
        dVdt.rightCols(14) = A() * V.rightCols(14);
    }
};

class ode_dVdt{
    double sigma;
    Vector3d u_t, u_t1;
public:
    ode_dVdt(Vector3d u_t, Vector3d u_t1, double sigma) : u_t(u_t), u_t1(u_t1), sigma(sigma){}

    void operator()(const state_type1 &V, state_type1 &dVdt, const double t){
        dVdt.setZero();
        Vector3d u = u_t + t / dt * (u_t1 - u_t);
        A.Update(V.col(0), u, sigma);
        B.Update(V.col(0), u, sigma);
        f.Update(V.col(0), u);

        double alpha = t / dt;
        double beta = 1 - alpha;
        MatrixXd Phi_A_xi(14,15);
        Phi_A_xi.col(0) = V.col(0);
        Phi_A_xi.rightCols(14).setIdentity();

        ode_dPhidt dPhidt(u_t, u_t1, sigma);
        runge_kutta_dopri5<state_type2, double, state_type2, double, vector_space_algebra> stepper2;
        integrate_adaptive(make_controlled( 1E-12 , 1E-12 , stepper2 ), dPhidt, Phi_A_xi, t, dt, dt/5.);

        dVdt.block<14, 1>(0, 0) = sigma * f();
        dVdt.block<14, 14>(0, 1) = A()* Phi_A_xi.block<14,14>(0,1);
        dVdt.block<14, 3>(0, 15) = Phi_A_xi.block<14,14>(0,1) * B() * alpha;
        dVdt.block<14, 3>(0, 18) = Phi_A_xi.block<14,14>(0,1) * B() * beta;
        dVdt.block<14, 1>(0, 21) = Phi_A_xi.block<14,14>(0,1) * f();
        dVdt.block<14, 1>(0, 22) = Phi_A_xi.block<14,14>(0,1) * (-A() * V.col(0) - B() * u);

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
        X.col(k).segment(7, 4) << 1, 0, 0, 0;
        X.col(k).segment(11, 3) = alpha1 * x_init.segment(11, 3) + alpha2 * x_final.segment(11, 3);

        U.col(k) = X(0, k) * -g_I;
    }

    cout << "Initialization finished." << endl;

//START SUCCESSIVE CONVEXIFICATION
    double sigma = sigma_guess;
    state_type1 V;

    for(int it = 1; it < iterations + 1; it++) {
        cout << "Iteration " << it << endl;
        cout << "Calculating new transition matrices." << endl;
        const clock_t begin_time = clock();

        for (int k = 0; k < K-1; k++) {
            V.setZero();
            V.col(0) = X.col(k);
            V.block<14,14>(0, 1).setIdentity();

            ode_dVdt dVdt(U.col(k), U.col(k+1), sigma);
            runge_kutta_dopri5<state_type1, double, state_type1, double, vector_space_algebra> stepper1;
            integrate_adaptive(make_controlled(1E-12 , 1E-12 , stepper1), dVdt, V, 0., dt, dt/5.);

        }
        cout << "Transition matrices calculated in " << double( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds." << endl;


    }

}