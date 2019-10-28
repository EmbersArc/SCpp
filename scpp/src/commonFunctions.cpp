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

} // namespace scpp