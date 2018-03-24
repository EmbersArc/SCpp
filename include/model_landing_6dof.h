#pragma once


#include <Eigen/Dense>
using namespace Eigen;


class model_landing_6dof {

    const Vector3d g_I = Vector3d(-1, 0, 0);
    const Vector3d J_B = Vector3d(1e-2, 1e-2, 1e-2);
    const Vector3d r_T_B = Vector3d(-1e-2, 0, 0);
    const double alpha_m = 0.02;

public:

    using StateVector   = Eigen::Matrix<double, 14,  1>;
    using ControlVector = Eigen::Matrix<double,  3,  1>;
    using StateMatrix   = Eigen::Matrix<double, 14, 14>;
    using ControlMatrix = Eigen::Matrix<double, 14,  3>;

    void initialize(int K, MatrixXd &X, MatrixXd &U);

    StateVector                ode(const StateVector &x, const ControlVector &u);
    StateMatrix     state_jacobian(const StateVector &x, const ControlVector &u);
    ControlMatrix control_jacobian(const StateVector &x, const ControlVector &u);
    

};