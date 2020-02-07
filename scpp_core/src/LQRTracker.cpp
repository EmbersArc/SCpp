#include "LQRTracker.hpp"

namespace scpp
{

LQRTracker::LQRTracker(Model::ptr_t model, const trajectory_data_t &td)
    : model(model), td(td)
{
    loadParameters();

    gains.resize(td.n_X());
    for (size_t k = 0; k < td.n_X(); k++)
    {
        Model::state_matrix_t A;
        Model::control_matrix_t B;

        if (not td.interpolatedInput() and k == td.n_X() - 2)
        {
            model->computeJacobians(td.X.at(k), td.U.at(k - 1), A, B);
        }
        else
        {
            model->computeJacobians(td.X.at(k), td.U.at(k), A, B);
        }

        ComputeLQR(Q, R, A, B, gains.at(k));
    }
}

void LQRTracker::loadParameters()
{
    ParameterServer param(model->getParameterFolder() + "/LQR.info");

    Model::state_vector_t q;
    Model::input_vector_t r;
    param.loadMatrix("state_weights", q);
    param.loadMatrix("input_weights", r);

    Q = q.asDiagonal();
    R = r.asDiagonal();
}

void LQRTracker::getInput(double t, const Model::state_vector_t &x, Model::input_vector_t &u) const
{
    t = std::clamp(t, 0., td.t);
    const Model::state_vector_t x_target = td.approxStateAtTime(t);
    const Model::input_vector_t u_target = td.inputAtTime(t);
    const Model::feedback_matrix_t K = interpolateGains(t);

    u = -K * (x - x_target) + u_target;
}

Model::feedback_matrix_t LQRTracker::interpolateGains(double t) const
{
    t = std::clamp(t, 0., td.t);

    const double dt = td.t / (td.n_X() - 1);
    double interpolate_value = std::fmod(t, dt) / dt;
    const size_t i = t / dt;

    const Model::feedback_matrix_t K0 = gains.at(i);
    const Model::feedback_matrix_t K1 = td.interpolatedInput() ? gains.at(std::min(gains.size() - 1, i + 1)) : gains.at(i);

    return K0 + interpolate_value * (K1 - K0);
}

} // namespace scpp