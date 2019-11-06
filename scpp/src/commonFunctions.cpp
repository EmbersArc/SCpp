#include "commonFunctions.hpp"

namespace scpp
{

Model::input_vector_t interpolatedInput(const Model::input_vector_v_t &U, double t,
                                        double total_time, bool first_order_hold)
{
    const size_t K = U.size();
    const double time_step = total_time / (K - 1);
    const size_t i = std::min(size_t(t / time_step), K - 2);
    const Model::input_vector_t u0 = U.at(i);
    const Model::input_vector_t u1 = first_order_hold ? U.at(i + 1) : u0;

    const double t_intermediate = std::fmod(t, time_step) / time_step;
    const Model::input_vector_t u = u0 + (u1 - u0) * t_intermediate;

    return u;
}

double expMovingAverage(double previousAverage, double period, double newValue)
{
    const double factor = 2. / (period + 1.);
    const double result = (newValue - previousAverage) * factor + previousAverage;
    return result;
}

std::vector<Eigen::Vector3d> getAccelerationRotatingFrame(const trajectory_data_t &td,
                                                          const Eigen::Vector3d offset)
{
    // calculate accelerations
    std::vector<Eigen::Vector3d> acc_passenger_b;
    const double dt = td.t / (td.n_X() - 1);

    for (size_t k = 0; k < td.n_X() - 1; k++)
    {
        const Eigen::Vector3d r_p_b = offset;
        const Eigen::Vector3d v0 = td.X.at(k).segment<3>(4);
        const Eigen::Vector3d v1 = td.X.at(k + 1).segment<3>(4);

        Eigen::Quaterniond q0;
        q0.w() = td.X.at(k)(7);
        q0.vec() << td.X.at(k).segment<3>(8);
        q0.normalize();
        Eigen::Quaterniond q1;
        q1.w() = td.X.at(k + 1)(7);
        q1.vec() << td.X.at(k + 1).segment<3>(8);
        q1.normalize();
        const Eigen::Quaterniond q = q0.slerp(0.5, q1);

        const Eigen::Vector3d w0 = td.X.at(k).segment<3>(11);
        const Eigen::Vector3d w1 = td.X.at(k + 1).segment<3>(11);
        const Eigen::Vector3d w = (w1 - w0) / 2;

        const Eigen::Vector3d vp0_i = v0 + w0.cross(q0 * r_p_b);
        const Eigen::Vector3d vp1_i = v1 + w1.cross(q1 * r_p_b);

        const Eigen::Vector3d dw = (w1 - w0) / dt;
        const Eigen::Vector3d a_i = (vp1_i - vp0_i) / dt;
        const Eigen::Vector3d g(0., 0., 0.);

        const Eigen::Vector3d a_b = a_i - w.cross(w.cross(r_p_b)) - dw.cross(r_p_b) + q.inverse() * g;
        acc_passenger_b.push_back(a_b);
    }
    return acc_passenger_b;
}

} // namespace scpp