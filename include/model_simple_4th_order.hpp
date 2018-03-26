#pragma once


#include <Eigen/Dense>
using namespace Eigen;


class model_simple_4th_order {
public:

    static constexpr size_t n_states = 4;
    static constexpr size_t n_inputs = 1;

    using StateVector   = Eigen::Matrix<double, n_states,        1>;
    using ControlVector = Eigen::Matrix<double, n_inputs,        1>;
    using StateMatrix   = Eigen::Matrix<double, n_states, n_states>;
    using ControlMatrix = Eigen::Matrix<double, n_states, n_inputs>;

    double total_time_guess() { return 3; }

    void initialize(int K, MatrixXd &X, MatrixXd &U);

    StateVector                ode(const StateVector &x, const ControlVector &u);
    StateMatrix     state_jacobian(const StateVector &x, const ControlVector &u);
    ControlMatrix control_jacobian(const StateVector &x, const ControlVector &u);
    

};