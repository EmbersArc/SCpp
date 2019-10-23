#include "optimizationProblem.hpp"
#include "MPCProblem.hpp"

namespace scpp
{

op::SecondOrderConeProgram buildMPCProblem(
    Model::ptr_t model,
    Model::state_vector_v_t &X,
    Model::input_vector_v_t &U,
    Model::state_vector_t &x_init,
    Model::state_vector_t &x_final,
    Model::state_vector_t &state_weights_intermediate,
    Model::state_vector_t &state_weights_terminal,
    Model::input_vector_t &input_weights,
    Model::state_matrix_t &A,
    Model::control_matrix_t &B,
    Model::state_vector_t &z,
    bool constant_dynamics,
    bool intermediate_cost_active)
{
    const size_t K = X.size();

    op::SecondOrderConeProgram socp;

    socp.createTensorVariable("X", {Model::state_dim, K});     // states
    socp.createTensorVariable("U", {Model::input_dim, K - 1}); // inputs
    socp.createTensorVariable("error_cost");                   // error minimization term
    socp.createTensorVariable("input_cost");                   // input minimization term

    // shortcuts to access solver variables and create parameters
    auto var = [&socp](const std::string &name, const std::vector<size_t> &indices = {}) { return socp.getVariable(name, indices); };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    // Initial state
    for (size_t i = 0; i < Model::state_dim; i++)
    {
        socp.addConstraint((-1.0) * var("X", {i, 0}) + param(x_init(i)) == 0.0);
    }

    for (size_t k = 0; k < K - 1; k++)
    {
        /**
         * Build linearized model equality constraint
         *    x(k+1) == A x(k) + B u(k) + z
         * -I x(k+1)  + A x(k) + B u(k) + z == 0
         * 
         */
        for (size_t i = 0; i < Model::state_dim; i++)
        {
            // -I * x(k+1)
            op::AffineExpression eq = (-1.0) * var("X", {i, k + 1});

            // A * x(k)
            for (size_t j = 0; j < Model::state_dim; j++)
                if (constant_dynamics)
                {
                    if (A(i, j) != 0)
                        eq = eq + A(i, j) * var("X", {j, k});
                }
                else
                {
                    eq = eq + param(A(i, j)) * var("X", {j, k});
                }

            // B * u(k)
            for (size_t j = 0; j < Model::input_dim; j++)
                if (constant_dynamics)
                {
                    if (B(i, j) != 0.)
                        eq = eq + B(i, j) * var("U", {j, k});
                }
                else
                {
                    eq = eq + param(B(i, j)) * var("U", {j, k});
                }

            // z
            if (constant_dynamics)
            {
                if (z(i) != 0.)
                    eq = eq + z(i);
            }
            else
            {
                eq = eq + param(z(i));
            }

            socp.addConstraint(eq == 0.0);
        }
    }

    /**
     * Build error cost
     * 
     */
    std::vector<op::AffineExpression> error_norm2_args;
    if (intermediate_cost_active)
    {
        for (size_t k = 1; k < K - 1; k++)
        {
            for (size_t i = 0; i < Model::state_dim; i++)
            {
                op::Parameter x_desired =
                    param_fn([&state_weights_intermediate, &x_final, i]() { return -1.0 * state_weights_intermediate(i) * x_final(i); });
                op::AffineTerm x_current =
                    param(state_weights_intermediate(i)) * var("X", {i, k});
                op::AffineExpression ex = x_desired + x_current;
                error_norm2_args.push_back(ex);
            }
        }
    }
    for (size_t i = 0; i < Model::state_dim; i++)
    {
        op::Parameter x_desired =
            param_fn([&state_weights_terminal, &x_final, i]() { return -1.0 * state_weights_terminal(i) * x_final(i); });
        op::AffineTerm x_current =
            param(state_weights_terminal(i)) * var("X", {i, K - 1});
        op::AffineExpression ex = x_desired + x_current;
        error_norm2_args.push_back(ex);
    }
    socp.addConstraint(op::norm2(error_norm2_args) <= (1.0) * var("error_cost"));
    socp.addMinimizationTerm(1.0 * var("error_cost"));

    /**
     * Build input cost
     * 
     */
    std::vector<op::AffineExpression> input_norm2_args;
    for (size_t k = 0; k < K - 1; k++)
    {
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            op::AffineExpression ex = param(input_weights(i)) * var("U", {i, k});
            input_norm2_args.push_back(ex);
        }
    }
    socp.addConstraint(op::norm2(input_norm2_args) <= (1.0) * var("input_cost"));
    socp.addMinimizationTerm(1.0 * var("input_cost"));

    model->addApplicationConstraints(socp, X, U);
    return socp;
}

} // namespace scpp
