#pragma once


#include <Eigen/Dense>
using namespace Eigen;

class model_landing_6dof {
public:
    void initialize(int K, MatrixXd &X, MatrixXd &U);
};