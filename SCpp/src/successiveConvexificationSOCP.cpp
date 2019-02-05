#include "successiveConvexificationSOCP.hpp"

namespace sc
{

op::SecondOrderConeProgram buildSCSOCP(
    Model &model,
    double &weight_trust_region_time,
    double &weight_trust_region_trajectory,
    double &weight_virtual_control,
    Eigen::MatrixXd &X,
    Eigen::MatrixXd &U,
    double &sigma,
    Model::state_matrix_v_t &A_bar,
    Model::control_matrix_v_t &B_bar,
    Model::control_matrix_v_t &C_bar,
    Model::state_vector_v_t &S_bar,
    Model::state_vector_v_t &z_bar)
{
    const size_t K = X.cols();

    op::SecondOrderConeProgram socp;

    socp.createTensorVariable("X", {Model::state_dim_, K});            // states
    socp.createTensorVariable("U", {Model::input_dim_, K});            // inputs
    socp.createTensorVariable("nu", {Model::state_dim_, K - 1});       // virtual control
    socp.createTensorVariable("nu_bound", {Model::state_dim_, K - 1}); // virtual control
    socp.createTensorVariable("norm1_nu");                             // virtual control norm upper bound
    socp.createTensorVariable("sigma");                                // total time
    socp.createTensorVariable("Delta_sigma");                          // squared change of sigma
    socp.createTensorVariable("Delta", {K});                           // squared change of the stacked [ x(k), u(k) ] vector
    socp.createTensorVariable("norm2_Delta");                          // 2-norm of the Delta(k) variables

    // shortcuts to access solver variables and create parameters
    auto var = [&](const string &name, const vector<size_t> &indices = {}) { return socp.get_variable(name, indices); };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    // auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    // Main objective: minimize total time
    socp.addMinimizationTerm(1.0 * var("sigma"));

    for (size_t k = 0; k < K - 1; k++)
    {
        /*
         * Build linearized model equality constraint
         *    x(k+1) == A x(k) + B u(k) + C u(k+1) + Sigma sigma + z + nu
         * -I x(k+1)  + A x(k) + B u(k) + C u(k+1) + Sigma sigma + z + nu == 0
         * 
         */
        for (size_t row_index = 0; row_index < Model::state_dim_; row_index++)
        {
            // -I * x(k+1)
            op::AffineExpression eq = (-1.0) * var("X", {row_index, k + 1});

            // A * x(k)
            for (size_t col_index = 0; col_index < Model::state_dim_; ++col_index)
                eq = eq + param(A_bar.at(k)(row_index, col_index)) * var("X", {col_index, k});

            // B * u(k)
            for (size_t col_index = 0; col_index < Model::input_dim_; ++col_index)
                eq = eq + param(B_bar.at(k)(row_index, col_index)) * var("U", {col_index, k});

            // C * u(k+1)
            for (size_t col_index = 0; col_index < Model::input_dim_; ++col_index)
                eq = eq + param(C_bar.at(k)(row_index, col_index)) * var("U", {col_index, k + 1});

            // Sigma sigma
            eq = eq + param(S_bar.at(k)(row_index, 0)) * var("sigma");

            // z
            eq = eq + param(z_bar.at(k)(row_index, 0));

            // nu
            eq = eq + (1.0) * var("nu", {row_index, k});

            socp.addConstraint(eq == 0.0);
        }
    }

    /*
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
            for (size_t i = 0; i < Model::state_dim_; i++)
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

    /* 
     * Build sigma trust region
     * 
     *  norm2( sigma0 - sigma ) <= Delta_sigma;
     * 
     */
    {
        vector<op::AffineExpression> norm2_args;
        norm2_args = {param(sigma) + (-1.0) * var("sigma")};

        socp.addConstraint(op::norm2(norm2_args) <= (1.0) * var("Delta_sigma"));

        // Minimize Delta_sigma
        socp.addMinimizationTerm(param(weight_trust_region_time) * var("Delta_sigma"));
    }

    for (size_t k = 0; k < K; k++)
    {
        /* 
         * Build state and input trust-region:
         * 
         *  norm2( [(x - x0)^T | (u - u0)^T]^T ) <= Delta;
         * 
         */

        vector<op::AffineExpression> norm2_args;
        for (size_t i = 0; i < Model::state_dim_; i++)
        {
            norm2_args.push_back(param(X(i, k)) + (-1.0) * var("X", {i, k}));
        }
        for (size_t i = 0; i < Model::input_dim_; i++)
        {
            norm2_args.push_back(param(U(i, k)) + (-1.0) * var("U", {i, k}));
        }
        socp.addConstraint(op::norm2(norm2_args) <= (1.0) * var("Delta", {k}));
    }

    /*
     * Build combined state/input trust region over all K:
     *  
     *  norm2([ Delta(1), Delta(2), ... , Delta(K) ]) <= norm2_Delta
     * 
     */
    {
        vector<op::AffineExpression> norm2_args;
        for (size_t k = 0; k < K; k++)
        {
            norm2_args.push_back((1.0) * var("Delta", {k}));
        }
        socp.addConstraint(op::norm2(norm2_args) <= (1.0) * var("norm2_Delta"));

        // Minimize norm2_Delta
        socp.addMinimizationTerm(param(weight_trust_region_trajectory) * var("norm2_Delta"));
    }

    // Total time must not be negative
    socp.addConstraint((1.0) * var("sigma") + (-0.01) >= (0.0));

    model.addApplicationConstraints(socp, X, U);
    return socp;
}

} // namespace sc
