#include "model_landing_6dof.h"

#include <iostream>
#include <cmath>
#include <ctime>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using namespace std;
using namespace boost::numeric::odeint;


using Model = model_landing_6dof;
model_landing_6dof model;

//trajectory points
const int K = 50;
double dt = 1 / double(K-1);

const int iterations = 15;

typedef Eigen::Matrix<double,14,23> state_type;

class ode_dVdt{
private:
    Vector3d u_t, u_t1;
    double sigma;

public:
    void Update(const Vector3d &u_t, const Vector3d &u_t1, const double &sigma){
        this->u_t = u_t;
        this->u_t1 = u_t1;
        this->sigma = sigma;
    }

    void operator()(const state_type &V, state_type &dVdt, const double t){

        const Model::StateVector &x = V.col(0);
        const Vector3d u = u_t + t / dt * (u_t1 - u_t);

        const double alpha = t / dt;
        const double beta = 1. - alpha;

        const Model::StateMatrix   A_bar  = sigma * model.state_jacobian(x, u);
        const Model::ControlMatrix B_bar  = sigma * model.control_jacobian(x, u);
        const Model::StateVector   f      =         model.ode(x, u);


        Model::StateMatrix Phi_A_xi = V.block<14, 14>(0, 1);
        Model::StateMatrix Phi_A_xi_inverse = Phi_A_xi.inverse();

        dVdt.block<14, 1>(0, 0) = sigma * f;
        dVdt.block<14, 14>(0, 1) = A_bar * Phi_A_xi;
        dVdt.block<14, 3>(0, 15) = Phi_A_xi_inverse * B_bar * alpha;
        dVdt.block<14, 3>(0, 18) = Phi_A_xi_inverse * B_bar * beta;
        dVdt.block<14, 1>(0, 21) = Phi_A_xi_inverse * f;
        dVdt.block<14, 1>(0, 22) = Phi_A_xi_inverse * (-A_bar * x - B_bar * u);

    }
};

int main() {

    MatrixXd X(14, K);
    MatrixXd U(3, K);


//START INITIALIZATION
    cout << "Starting initialization." << endl;

    model.initialize(K, X, U);

    cout << "Initialization finished." << endl;

//START SUCCESSIVE CONVEXIFICATION
    
    const double sigma_guess = 3.;// TODO from model
    double sigma = sigma_guess;

    state_type V;

    // TODO array<> instead of vector<>
    vector<Model::StateMatrix> A_bar(K);
    vector<Model::ControlMatrix> B_bar(K);
    vector<Model::ControlMatrix> C_bar(K);
    vector<Model::StateVector> Sigma_bar(K);
    vector<Model::StateVector> z_bar(K);

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