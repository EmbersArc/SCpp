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
                                                          const Eigen::Vector3d offset,
                                                          const double g)
{
    Model::state_vector_v_t X = td.X;
    X.push_back(X.back());

    // calculate accelerations
    std::vector<Eigen::Vector3d> acc_passenger_b;
    const double dt = td.t / (X.size() - 1);

    for (size_t k = 0; k < X.size() - 1; k++)
    {
        const Eigen::Vector3d r_p_b = offset;
        const Eigen::Vector3d v0 = X.at(k).segment<3>(4);
        const Eigen::Vector3d v1 = X.at(k + 1).segment<3>(4);

        // Eigen::Quaterniond q0;
        // q0.w() = X.at(k)(7);
        // q0.vec() << X.at(k).segment<3>(8);
        // q0.normalize();
        // Eigen::Quaterniond q1;
        // q1.w() = X.at(k + 1)(7);
        // q1.vec() << X.at(k + 1).segment<3>(8);
        // q1.normalize();
        // const Eigen::Quaterniond q = q0.slerp(0.5, q1);

        const Eigen::Vector3d w0 = X.at(k).segment<3>(11);
        const Eigen::Vector3d w1 = X.at(k + 1).segment<3>(11);
        const Eigen::Vector3d w_b = (w1 - w0) / 2;

        const Eigen::Vector3d dw_b = (w1 - w0) / dt;
        const Eigen::Vector3d dv_i = (v1 - v0) / dt;

        const Eigen::Vector3d a_coriolis(0., 0., 0.); // r_p_b is constant
        const Eigen::Vector3d a_centrifugal = -w_b.cross(w_b.cross(r_p_b));
        const Eigen::Vector3d a_euler = -dw_b.cross(r_p_b);
        const Eigen::Vector3d a_imp = dv_i + Eigen::Vector3d(0., 0., g);

        acc_passenger_b.push_back(a_imp + a_centrifugal + a_coriolis + a_euler);
    }
    return acc_passenger_b;
}

std::string getTimeString(){
    using sc = std::chrono::system_clock ;
    std::time_t t = sc::to_time_t(sc::now());
    char buf[20];
    std::strftime(buf, 20, "%Y_%m_%d_%H_%M_%S", std::localtime(&t));
    return std::string(buf);
}

} // namespace scpp