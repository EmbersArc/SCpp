#pragma once

#include "LQR.hpp"
#include "parameterServer.hpp"

namespace scpp
{

class LQRTracker
{
private:
    std::vector<Model::feedback_matrix_t> gains;

    Model::state_matrix_t Q;
    Model::input_matrix_t R;

    Model::ptr_t model;
    trajectory_data_t td;

    void loadParameters();

    Model::feedback_matrix_t interpolateGains(double t) const;

public:
    LQRTracker(Model::ptr_t model, const trajectory_data_t &td);
    void getInput(double t, const Model::state_vector_t &x, Model::input_vector_t &u) const;
};

} // namespace scpp