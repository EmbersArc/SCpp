#pragma once

#include <cstddef>

namespace scpp
{

template <class Model>
struct TrajectoryData
{
    typename Model::state_vector_v_t X;
    typename Model::input_vector_v_t U;
    double t;

    void initialize(size_t K, bool interpolate_input);

    bool interpolatedInput() const;

    size_t n_X() const;
    size_t n_U() const;
};

template <class Model>
void TrajectoryData<Model>::initialize(size_t K, bool interpolate_input)
{
    X.resize(K);
    U.resize(interpolate_input ? K : K - 1);
    t = 0.;
}

template <class Model>
bool TrajectoryData<Model>::interpolatedInput() const
{
    return U.size() == X.size();
}

template <class Model>
size_t TrajectoryData<Model>::n_X() const
{
    return X.size();
}

template <class Model>
size_t TrajectoryData<Model>::n_U() const
{
    return U.size();
}

} // namespace scpp