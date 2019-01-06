#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

#include "active_model.hpp"

class ODE{

private:
    Model model;
    Model::ControlVector U0;
    Model::ControlVector U1;
    double dt;

public:
    explicit ODE(double dt):dt(dt){};

    void operator()(const Model::StateVector &X, Model::StateVector &dXdt, double t);
    void update_input(const Model::ControlVector &U0, const Model::ControlVector &U1);

};

void check_trajectory(
        Eigen::Matrix<double, Model::n_states, K> &X_in,
        Eigen::Matrix<double, Model::n_inputs, K> &U,
        double sigma
);