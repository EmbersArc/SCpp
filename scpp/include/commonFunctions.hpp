#pragma once

#include "activeModel.hpp"

namespace scpp
{

Model::input_vector_t interpolatedInput(const Model::input_vector_v_t &U, double t,
                                        double total_time, bool first_order_hold);

double expMovingAverage(double previousAverage, double period, double newValue);

} // namespace scpp