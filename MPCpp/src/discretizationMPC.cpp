#include <boost/numeric/odeint.hpp>

#include "eigenIntegration.hpp"
#include <eigen3/unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h>

#include "discretizationMPC.hpp"

void eulerLinearDiscretization(Model &model,
                               double &ts,
                               const Model::state_vector_t &x_eq,
                               const Model::input_vector_t &u_eq,
                               Model::state_matrix_t &A,
                               Model::control_matrix_t &B)
{
    Model::state_matrix_t A_c;
    Model::control_matrix_t B_c;
    model.computeJacobians(x_eq, u_eq, A_c, B_c);

    A = Model::state_matrix_t::Identity() + ts * A_c;
    B = ts * B_c;
}

class DiscretizationODE
{
private:
    Model::state_matrix_t A_c;
    Model::control_matrix_t B_c;
    double ts;

public:
    DiscretizationODE(
        Model::state_matrix_t &A_c,
        Model::control_matrix_t &B_c,
        double ts) : A_c(A_c), B_c(B_c), ts(ts)
    {
    }

    void operator()(const Model::control_matrix_t &B, Model::control_matrix_t &Bdt, const double t) const
    {
        Bdt = (A_c * (ts - t)).exp() * B_c;
    }
};

void exactLinearDiscretization(Model &model,
                               double &ts,
                               const Model::state_vector_t &x_eq,
                               const Model::input_vector_t &u_eq,
                               Model::state_matrix_t &A,
                               Model::control_matrix_t &B)
{
    Model::state_matrix_t A_c;
    Model::control_matrix_t B_c;
    model.computeJacobians(x_eq, u_eq, A_c, B_c);

    using namespace boost::numeric::odeint;
    runge_kutta4<Model::control_matrix_t, double, Model::control_matrix_t, double, vector_space_algebra> stepper;
    DiscretizationODE discretizationODE(A_c, B_c, ts);

    A = (ts * A_c).exp();
    B.setZero();
    integrate_adaptive(stepper, discretizationODE, B, 0., ts, ts / 10.);
}
