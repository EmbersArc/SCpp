#pragma once

#include "activeModel.hpp"

namespace scpp
{

Model::input_vector_t interpolatedInput(const Model::input_vector_v_t &U, double t,
                                        double total_time, bool first_order_hold);

double expMovingAverage(double previousAverage, double period, double newValue);

std::vector<Eigen::Vector3d> getAccelerationRotatingFrame(const trajectory_data_t &td,
                                                          const Eigen::Vector3d offset = Eigen::Vector3d::Zero(),
                                                          const double g = 0.);

std::string getTimeString(){
    using sc = std::chrono::system_clock ;
    std::time_t t = sc::to_time_t(sc::now());
    char buf[20];
    std::strftime(buf, 20, "%Y_%m_%d_%H_%M_%S", std::localtime(&t));
    return std::string(buf);
}

} // namespace scpp