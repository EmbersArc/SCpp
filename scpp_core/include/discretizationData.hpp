#pragma once

#include <cstddef>

namespace scpp
{

template<class Model>
struct DiscretizationData
{
    typename Model::state_matrix_v_t A;
    typename Model::control_matrix_v_t B;
    typename Model::control_matrix_v_t C;
    typename Model::state_vector_v_t s;
    typename Model::state_vector_v_t z;

    void initialize(size_t K, bool interpolate_input, bool free_final_time);

    bool interpolatedInput() const;

    bool variableTime() const;

    size_t n_X() const;
    size_t n_U() const;
};

template<class Model>
void DiscretizationData<Model>::initialize(size_t K, bool interpolate_input, bool free_final_time)
{
    A.resize(K - 1);

    B.resize(K - 1);

    if (interpolate_input)
    {
        C.resize(K - 1);
    }
    else
    {
        C.clear();
    }

    if (free_final_time)
    {
        s.resize(K - 1);
    }
    else
    {
        s.clear();
    }

    z.resize(K - 1);
}

template<class Model>
bool DiscretizationData<Model>::interpolatedInput() const
{
    return not C.empty();
}

template<class Model>
bool DiscretizationData<Model>::variableTime() const
{
    return not s.empty();
}

template<class Model>
size_t DiscretizationData<Model>::n_X() const
{
    return A.size();
}

template<class Model>
size_t DiscretizationData<Model>::n_U() const
{
    return B.size();
}

} // namespace scpp