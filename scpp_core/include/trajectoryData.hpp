#pragma once

#include "activeModel.hpp"

struct TrajectoryData
{
    Model::state_matrix_v_t A;
    Model::control_matrix_v_t B;
    Model::control_matrix_v_t C;
    Model::state_vector_v_t s;
    Model::state_vector_v_t z;

    void initialize(size_t K, bool interpolate_input, bool free_final_time);

    bool interpolatedInput() const;

    bool variableTime() const;
};
