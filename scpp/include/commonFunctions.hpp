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

std::string getTimeString();

template <typename T>
std::vector<T> reduce_vector(const std::vector<T> &v, size_t steps)
{
    const size_t size = v.size();

    std::vector<T> new_vector;

    for (size_t i = 0; i < steps; i++)
    {
        const size_t index = size_t(size / steps * i);
        new_vector.push_back(v.at(index));
    }
    return new_vector;
}

} // namespace scpp