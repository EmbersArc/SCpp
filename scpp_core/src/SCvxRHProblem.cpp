#include "SCvxRHProblem.hpp"
#include "SCvxProblem.hpp"

namespace scpp
{

op::SecondOrderConeProgram buildSCvxRHProblem(
    double &trust_region,
    double &weight_virtual_control,
    Model::state_vector_t &state_weights,
    Model::input_vector_t &input_weights,
    Model::state_vector_t &x_init,
    Model::state_vector_t &x_final,
    trajectory_data_t &td,
    discretization_data_t &dd)
{
    op::SecondOrderConeProgram socp = buildSCvxProblem(trust_region, weight_virtual_control, td, dd);

    // shortcuts to access solver variables and create parameters
    auto var = [&socp](const std::string &name, const std::vector<size_t> &indices = {}) { return socp.getVariable(name, indices); };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    // Initial state
    for (size_t i = 0; i < Model::state_dim; i++)
    {
        socp.addConstraint((-1.0) * var("X", {i, 0}) + param(x_init(i)) == 0.0);
    }

    /**
     * Build error cost
     * 
     */
    socp.createTensorVariable("error_cost"); // error minimization term
    std::vector<op::AffineExpression> error_norm2_args;
    for (size_t i = 0; i < Model::state_dim; i++)
    {
        op::Parameter x_desired =
            param_fn([&state_weights, &x_final, i]() { return -1.0 * state_weights(i) * x_final(i); });
        op::AffineTerm x_current =
            param(state_weights(i)) * var("X", {i, td.n_X() - 1});
        op::AffineExpression ex = x_desired + x_current;
        error_norm2_args.push_back(ex);
    }
    socp.addConstraint(op::norm2(error_norm2_args) <= (1.0) * var("error_cost"));
    socp.addMinimizationTerm(1.0 * var("error_cost"));

    /**
     * Build input cost
     * 
     */
    socp.createTensorVariable("input_cost"); // input minimization term
    std::vector<op::AffineExpression> input_norm2_args;
    for (size_t k = 0; k < td.n_U(); k++)
    {
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            op::AffineExpression ex = param(input_weights(i)) * var("U", {i, k});
            input_norm2_args.push_back(ex);
        }
    }
    socp.addConstraint(op::norm2(input_norm2_args) <= (1.0) * var("input_cost"));
    socp.addMinimizationTerm(1.0 * var("input_cost"));

    return socp;
}

} // namespace scpp
