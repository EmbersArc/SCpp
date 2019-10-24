#pragma once

#include "activeModel.hpp"

struct DiscretizationData
{
    Model::state_matrix_v_t A;
    Model::control_matrix_v_t B;
    Model::control_matrix_v_t C;
    Model::state_vector_v_t s;
    Model::state_vector_v_t z;

    bool interpolatedInput() const
    {
        return not C.empty();
    }

    bool variableTime() const
    {
        return not s.empty();
    }

    void initialize(size_t K, bool interpolate_input, bool free_final_time)
    {
        A.resize(K - 1);
        B.resize(K - 1);
        if (interpolate_input)
        {
            C.resize(K - 1);
        }
        if (free_final_time)
        {
            s.resize(K - 1);
        }
        z.resize(K - 1);
    }
};