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

    typename Model::input_vector_t inputAtTime(double t) const;
    typename Model::state_vector_t approxStateAtTime(double t) const;

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
typename Model::input_vector_t TrajectoryData<Model>::inputAtTime(double t) const
{
    t = std::clamp(t, 0., this->t);

    if (t == this->t)
    {
        return U.back();
    }

    const double dt = this->t / (n_X() - 1);
    double interpolate_value = std::fmod(t, dt) / dt;
    const size_t i = t / dt;

    const typename Model::input_vector_t u0 = U.at(i);
    const typename Model::input_vector_t u1 = interpolatedInput() ? U.at(i + 1) : U.at(i);

    return u0 + interpolate_value * (u1 - u0);
}

template <class Model>
typename Model::state_vector_t TrajectoryData<Model>::approxStateAtTime(double t) const
{
    t = std::clamp(t, 0., this->t);

    if (t == this->t)
    {
        return X.back();
    }

    const double dt = this->t / (n_X() - 1);
    double interpolate_value = std::fmod(t, dt) / dt;
    const size_t i = t / dt;

    const typename Model::state_vector_t x0 = X.at(i);
    const typename Model::state_vector_t x1 = X.at(i + 1);

    return x0 + interpolate_value * (x1 - x0);
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