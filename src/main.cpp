#include "params.h"
#include "model_landing_6dof.h"

AMatrix A;

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

        alpha = t / dt;
        beta = 1. - alpha;

        Phi_A_xi = V.block<14, 14>(0, 1).inverse();

        dVdt.block<14, 1>(0, 0) = sigma * A.get_f();
        dVdt.block<14, 14>(0, 1) = A.get_A()* V.block<14, 14>(0, 1);
        dVdt.block<14, 3>(0, 15) = Phi_A_xi * A.get_B() * alpha;
        dVdt.block<14, 3>(0, 18) = Phi_A_xi * A.get_B() * beta;
        dVdt.block<14, 1>(0, 21) = Phi_A_xi * A.get_f();
        dVdt.block<14, 1>(0, 22) = Phi_A_xi * (-A.get_A() * V.col(0) - A.get_B() * u);

    }
};

int main() {

    model_landing_6dof model;

    MatrixXd X(14, K);
    MatrixXd U(3, K);


//START INITIALIZATION
    cout << "Starting initialization." << endl;

    model.initialize(K, X, U);

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

            A_bar[k] = V.block<14,14>(0, 1);
            B_bar[k] = V.block<14,14>(0, 1) * V.block<14,3>(0, 15);
            C_bar[k] = V.block<14,14>(0, 1) * V.block<14,3>(0, 18);
            Sigma_bar[k] = V.block<14,14>(0, 1) * V.block<14,1>(0, 21);
            z_bar[k] = V.block<14,14>(0, 1) * V.block<14,1>(0, 22);

            // debug print for refactoring, remove later
            cout << A_bar[k] << endl;
            cout << B_bar[k] << endl;
            cout << C_bar[k] << endl;
            cout << Sigma_bar[k] << endl;
            cout << z_bar[k] << endl;
        }
        //cout << "Transition matrices calculated in " << double( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds." << endl;

        // TODO: Solve problem.


    }

}