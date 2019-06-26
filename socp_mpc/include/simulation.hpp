#include <activeModel.hpp>

namespace sim
{

class ODE
{
public:
    ODE(Model &model, double dt, const Model::input_vector_t &u0, const Model::input_vector_t &u1);
    void operator()(const Model::state_vector_t &f, Model::state_vector_t &dfdt, const double t);

private:
    Model &model;
    Model::input_vector_t u0, u1;
    double dt;
};

void simulate(Model &model, double dt,
              const Model::state_vector_t &x0,
              const Model::input_vector_t &u0,
              const Model::input_vector_t &u1,
              Model::state_vector_t &x1);

} // namespace sim