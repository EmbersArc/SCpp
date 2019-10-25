#include "optimizationProblem.hpp"
#include "SCvxProblem.hpp"

namespace scpp
{

op::SecondOrderConeProgram buildSCvxProblem(
    double &trust_region,
    double &weight_virtual_control,
    Model::state_vector_v_t &X,
    Model::input_vector_v_t &U,
    DiscretizationData &dd)
{
    op::SecondOrderConeProgram socp;

    // shortcuts to access solver variables and create parameters
    auto var = [&socp](const std::string &name, const std::vector<size_t> &indices = {}) { return socp.getVariable(name, indices); };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    // auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    socp.createTensorVariable("X", {Model::state_dim, X.size()});            // states
    socp.createTensorVariable("U", {Model::input_dim, U.size()});            // inputs
    socp.createTensorVariable("nu", {Model::state_dim, X.size() - 1});       // virtual control
    socp.createTensorVariable("nu_bound", {Model::state_dim, X.size() - 1}); // virtual control lower/upper bound
    socp.createTensorVariable("norm1_nu");                                   // virtual control norm

    for (size_t k = 0; k < X.size() - 1; k++)
    {
        /**
         * Build linearized model equality constraint
         *    x(k+1) == A x(k) + B u(k) + C u(k+1) + z + nu
         *   -x(k+1)  + A x(k) + B u(k) + C u(k+1) + z + nu == 0
         * 
         */
        for (size_t i = 0; i < Model::state_dim; i++)
        {
            // - x(k+1)
            op::AffineExpression eq = (-1.0) * var("X", {i, k + 1});

            // A * x(k)
            for (size_t j = 0; j < Model::state_dim; j++)
                eq = eq + param(dd.A.at(k)(i, j)) * var("X", {j, k});

            // B * u(k)
            for (size_t j = 0; j < Model::input_dim; j++)
                eq = eq + param(dd.B.at(k)(i, j)) * var("U", {j, k});

            if (dd.interpolatedInput())
            {
                // C * u(k+1)
                for (size_t j = 0; j < Model::input_dim; j++)
                    eq = eq + param(dd.C.at(k)(i, j)) * var("U", {j, k + 1});
            }

            // z
            eq = eq + param(dd.z.at(k)(i, 0));

            // nu
            eq = eq + (1.0) * var("nu", {i, k});

            socp.addConstraint(eq == 0.0);
        }
    }

    /**
     * Build virtual control norm
     *
     * minimize (weight_virtual_control * norm1_nu)
     * s.t. sum(nu_bound) <= norm1_nu
     *      -nu_bound <= nu <= nu_bound
     *
     */
    {
        op::AffineExpression bound_sum;
        for (size_t k = 0; k < X.size() - 1; k++)
        {
            for (size_t i = 0; i < Model::state_dim; i++)
            {
                // -nu_bound <= nu
                socp.addConstraint((1.0) * var("nu_bound", {i, k}) + (1.0) * var("nu", {i, k}) >= (0.0));
                //  nu <= nu_bound
                socp.addConstraint((1.0) * var("nu_bound", {i, k}) + (-1.0) * var("nu", {i, k}) >= (0.0));

                // sum(-nu_bound)
                bound_sum = bound_sum + (-1.0) * var("nu_bound", {i, k});
            }
        }
        // sum(nu_bound) <= norm1_nu
        socp.addConstraint((1.0) * var("norm1_nu") + bound_sum >= (0.0));

        // Minimize the virtual control
        socp.addMinimizationTerm(param(weight_virtual_control) * var("norm1_nu"));
    }

    for (size_t k = 0; k < U.size(); k++)
    {
        /**
         * Build input trust region:
         *     norm2(u - u0)  <=  trust_region
         *
         */

        std::vector<op::AffineExpression> norm2_args;
        for (size_t i = 0; i < Model::input_dim; i++)
        {
            norm2_args.push_back(param(U[k](i)) + (-1.0) * var("U", {i, k}));
        }
        socp.addConstraint(op::norm2(norm2_args) <= param(trust_region));
    }

    return socp;
}

} // namespace scpp
