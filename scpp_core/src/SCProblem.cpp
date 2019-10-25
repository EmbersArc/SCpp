#include "optimizationProblem.hpp"
#include "SCProblem.hpp"

namespace scpp
{

op::SecondOrderConeProgram buildSCProblem(
    double &weight_time,
    double &weight_trust_region_time,
    double &weight_trust_region_trajectory,
    double &weight_virtual_control,
    Model::state_vector_v_t &X,
    Model::input_vector_v_t &U,
    double &sigma,
    DiscretizationData &dd)
{
    const size_t K = X.size();

    op::SecondOrderConeProgram socp;

    // shortcuts to access solver variables and create parameters
    auto var = [&socp](const std::string &name, const std::vector<size_t> &indices = {}) { return socp.getVariable(name, indices); };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    // auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    socp.createTensorVariable("X", {Model::state_dim, X.size()});            // states
    socp.createTensorVariable("U", {Model::input_dim, U.size()});            // inputs
    socp.createTensorVariable("nu", {Model::state_dim, X.size() - 1});       // virtual control
    socp.createTensorVariable("nu_bound", {Model::state_dim, X.size() - 1}); // virtual control
    socp.createTensorVariable("norm1_nu");                                   // virtual control norm upper bound
    socp.createTensorVariable("Delta", {K});                                 // squared change of the stacked [ x(k), u(k) ] vector
    socp.createTensorVariable("norm2_Delta");                                // 2-norm of the Delta(k) variables
    if (dd.variableTime())
    {
        socp.createTensorVariable("sigma");       // total time
        socp.createTensorVariable("Delta_sigma"); // squared change of sigma
        // minimize total time
        socp.addMinimizationTerm(param(weight_time) * var("sigma"));
        // Total time must not be negative
        socp.addConstraint((1.0) * var("sigma") + (-0.01) >= (0.0));
    }

    for (size_t k = 0; k < K - 1; k++)
    {
        /**
         * Build linearized model equality constraint
         *    x(k+1) == A x(k) + B u(k) + C u(k+1) + Sigma sigma + z + nu
         *   -x(k+1)  + A x(k) + B u(k) + C u(k+1) + Sigma sigma + z + nu == 0
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

            if (dd.variableTime())
            {
                // Sigma sigma
                eq = eq + param(dd.s.at(k)(i, 0)) * var("sigma");
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
        for (size_t k = 0; k < K - 1; k++)
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
        // sum(-nu_bound) <= norm1_nu
        socp.addConstraint((1.0) * var("norm1_nu") + bound_sum >= (0.0));

        // Minimize the virtual control
        socp.addMinimizationTerm(param(weight_virtual_control) * var("norm1_nu"));
    }

    if (dd.variableTime())
    {
        /**
        *  Build sigma trust region
        * (sigma - sigma0) * (sigma - sigma0) <= Delta_sigma
        *          is equivalent to
        * norm2(
        *        0.5 - 0.5 * Delta_sigma
        *        sigma0 - sigma
        *      )
        *      <= 0.5 + 0.5 * Delta_sigma;
        */
        {
            std::vector<op::AffineExpression> norm2_args;
            norm2_args = {(0.5) + (-0.5) * var("Delta_sigma"),
                          param(sigma) + (-1.0) * var("sigma")};

            socp.addConstraint(op::norm2(norm2_args) <= (0.5) + (0.5) * var("Delta_sigma"));

            // Minimize Delta_sigma
            socp.addMinimizationTerm(param(weight_trust_region_time) * var("Delta_sigma"));
        }
    }

    for (size_t k = 0; k < K; k++)
    {
        /**
         * Build state and input trust-region:
         *     (x - x0)^T * (x - x0)  +  (u - u0)^T * (u - u0)  <=  Delta
         *
         * The constraint is equivalent to the SOCP form:
         *
         * norm2(
         *        0.5 - 0.5 * Delta
         *        [(x - x0)^T | (u - u0)^T]^T
         *      )
         *     <= 0.5 + 0.5 * Delta;
         *
         */

        std::vector<op::AffineExpression> norm2_args;

        norm2_args.push_back((0.5) + (-0.5) * var("Delta", {k}));

        for (size_t i = 0; i < Model::state_dim; i++)
        {
            norm2_args.push_back(param(X[k](i)) + (-1.0) * var("X", {i, k}));
        }
        if (not(dd.interpolatedInput() and k == K - 1))
        {
            for (size_t i = 0; i < Model::input_dim; i++)
            {
                norm2_args.push_back(param(U[k](i)) + (-1.0) * var("U", {i, k}));
            }
        }
        socp.addConstraint(op::norm2(norm2_args) <= (0.5) + (0.5) * var("Delta", {k}));
    }

    /**
     * Build combined state/input trust region over all K:
     *  
     *  norm2([ Delta(1), Delta(2), ... , Delta(K) ]) <= norm2_Delta
     * 
     */
    {
        std::vector<op::AffineExpression> norm2_args;
        for (size_t k = 0; k < K; k++)
        {
            norm2_args.push_back((1.0) * var("Delta", {k}));
        }
        socp.addConstraint(op::norm2(norm2_args) <= (1.0) * var("norm2_Delta"));

        // Minimize norm2_Delta
        socp.addMinimizationTerm(param(weight_trust_region_trajectory) * var("norm2_Delta"));
    }

    return socp;
}

} // namespace scpp
