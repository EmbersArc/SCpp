#ifndef OPTIMALLANDING_CHECK_JACOBIAN_H
#define OPTIMALLANDING_CHECK_JACOBIAN_H

#include "model_landing_6dof.h"

#include <iostream>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>


using Model = model_landing_6dof;
using TestStateVector = Eigen::Matrix<double, 2 * Model::n_states, 1>;  // stacked state vector

class TestODE {
    Model model;
    Vector3d U0;
    Vector3d U1;

public:
    TestODE();
    void operator()(const TestStateVector &X, TestStateVector &dXdt, double t);
};

bool test_model(Model::StateVector V);

#endif //OPTIMALLANDING_CHECK_JACOBIAN_H
