#include <iostream>

#include "successiveConvexificationSOCP.hpp"

optimization_problem::SecondOrderConeProgram build_successive_convexification_SOCP(
    Model &model,
    double &weight_trust_region_sigma,
    double &weight_trust_region_xu,
    double &weight_virtual_control,
    Eigen::Matrix<double, Model::n_states, K> &X,
    Eigen::Matrix<double, Model::n_inputs, K> &U,
    double &sigma,
    array<Model::StateMatrix, (K - 1)> &A_bar,
    array<Model::ControlMatrix, (K - 1)> &B_bar,
    array<Model::ControlMatrix, (K - 1)> &C_bar,
    array<Model::StateVector, (K - 1)> &Sigma_bar,
    array<Model::StateVector, (K - 1)> &z_bar)
{

    const size_t n_states = Model::n_states;
    const size_t n_inputs = Model::n_inputs;

    optimization_problem::SecondOrderConeProgram socp;

    socp.create_tensor_variable("X", {n_states, K});            // states
    socp.create_tensor_variable("U", {n_inputs, K});            // inputs
    socp.create_tensor_variable("nu", {n_states, K - 1});       // virtual control
    socp.create_tensor_variable("nu_bound", {n_states, K - 1}); // virtual control
    socp.create_tensor_variable("norm1_nu", {});                // virtual control norm upper bound
    socp.create_tensor_variable("sigma", {});                   // total time
    socp.create_tensor_variable("Delta_sigma", {});             // squared change of sigma
    socp.create_tensor_variable("Delta", {K});                  // squared change of the stacked [ x(k), u(k) ] vector
    socp.create_tensor_variable("norm2_Delta", {});             // 2-norm of the Delta(k) variables

    // shortcuts to access solver variables and create parameters
    auto var = [&](const string &name, const vector<size_t> &indices) { return socp.get_variable(name, indices); };
    auto param = [](double &param_value) { return optimization_problem::Parameter(&param_value); };
    auto param_fn = [](std::function<double()> callback) { return optimization_problem::Parameter(callback); };

    // Main objective: minimize total time
    socp.add_minimization_term(1.0 * var("sigma", {}));

    for (size_t k = 0; k < K - 1; k++)
    {

        // Build linearized model equality constraint
        //    x(k+1) == A x(k) + B u(k) + C u(k+1) + Sigma sigma + z + nu
        // -I x(k+1)  + A x(k) + B u(k) + C u(k+1) + Sigma sigma + z + nu == 0
        for (size_t row_index = 0; row_index < n_states; ++row_index)
        {

            // -I * x(k+1)
            optimization_problem::AffineExpression eq = (-1.0) * var("X", {row_index, k + 1});

            // A * x(k)
            for (size_t col_index = 0; col_index < n_states; ++col_index)
                eq = eq + param(A_bar.at(k)(row_index, col_index)) * var("X", {col_index, k});

            // B * u(k)
            for (size_t col_index = 0; col_index < n_inputs; ++col_index)
                eq = eq + param(B_bar.at(k)(row_index, col_index)) * var("U", {col_index, k});

            // C * u(k+1)
            for (size_t col_index = 0; col_index < n_inputs; ++col_index)
                eq = eq + param(C_bar.at(k)(row_index, col_index)) * var("U", {col_index, k + 1});

            // Sigma sigma
            eq = eq + param(Sigma_bar.at(k)(row_index, 0)) * var("sigma", {});

            // z
            eq = eq + param(z_bar.at(k)(row_index, 0));

            // nu
            eq = eq + (1.0) * var("nu", {row_index, k});

            socp.add_constraint(eq == 0.0);
        }
    }

    // Build virtual control norm
    // minimize (weight_virtual_control * norm1_nu)
    // s.t. sum(nu_bound) <= norm1_nu
    //      -nu_bound <= nu <= nu_bound
    {
        optimization_problem::AffineExpression bound_sum;
        for (size_t k = 0; k < K - 1; k++)
        {
            for (size_t row_index = 0; row_index < n_states; ++row_index)
            {
                // -nu_bound <= nu
                socp.add_constraint((1.0) * var("nu_bound", {row_index, k}) + (1.0) * var("nu", {row_index, k}) >= (0.0));
                // nu <= nu_bound
                socp.add_constraint((1.0) * var("nu_bound", {row_index, k}) + (-1.0) * var("nu", {row_index, k}) >= (0.0));

                // sum(-nu_bound)
                bound_sum = bound_sum + (-1.0) * var("nu_bound", {row_index, k});
            }
        }
        // sum(-nu_bound) <= norm1_nu
        socp.add_constraint((1.0) * var("norm1_nu", {}) + bound_sum >= (0.0));

        // Minimize the virtual control
        socp.add_minimization_term(param(weight_virtual_control) * var("norm1_nu", {}));
    }

    // Build sigma trust region
    // (sigma-sigma0) * (sigma-sigma0) <= Delta_sigma
    //          is equivalent to
    //   norm2([
    //       ((-sigma0)*sigma   +(-0.5)*Delta_sigma  +(0.5+0.5*sigma0*sigma0)  ),
    //       sigma
    //   ])
    //   <= (  +(sigma0)*sigma  +(0.5)*Delta_sigma    +(0.5-0.5*sigma0*sigma0)   )
    {
        // Formulas involving sigma, from the above comment
        auto sigma_fn1 = param_fn([&sigma]() { return -sigma; });
        auto sigma_fn2 = param_fn([&sigma]() { return (0.5 + 0.5 * sigma * sigma); });
        auto sigma_fn3 = param_fn([&sigma]() { return sigma; });
        auto sigma_fn4 = param_fn([&sigma]() { return (0.5 - 0.5 * sigma * sigma); });

        socp.add_constraint(
            optimization_problem::norm2({sigma_fn1 * var("sigma", {}) + (-0.5) * var("Delta_sigma", {}) + sigma_fn2,
                                         (1.0) * var("sigma", {})}) <= sigma_fn3 * var("sigma", {}) + (0.5) * var("Delta_sigma", {}) + sigma_fn4);

        // Minimize Delta_sigma
        socp.add_minimization_term(param(weight_trust_region_sigma) * var("Delta_sigma", {}));
    }

    for (size_t k = 0; k < K; k++)
    {
        /* 
         * Build state and input trust-region:
         *     (x - x0)^T * (x - x0)  +  (u - u0)^T * (u - u0)  <=  Delta
         * the index k is omitted, but applies to all terms in the constraint.
         * The constraint is equivalent to the SOCP form:
         * 
         * norm2(
         *         ( (-x0^T)*x  +(-u0^T)*u  +(-0.5)*Delta  +(0.5 + 0.5*x0^T*x0 + 0.5*u0^T*u0)),
         *         (I)*x, 
         *         (I)*u
         * )
         *     <= (  ( x0^T)*x  +( u0^T)*u  +( 0.5)*Delta  +(0.5 - 0.5*x0^T*x0 - 0.5*u0^T*u0));
         * 
         */

        vector<optimization_problem::AffineExpression> norm2_args;

        { // (-x0^T)*x  +(-u0^T)*u  +(-0.5)*Delta  +(0.5 + 0.5*x0^T*x0 + 0.5*u0^T*u0)
            optimization_problem::AffineExpression norm2_first_arg;

            // (-x0^T)*x
            for (size_t i = 0; i < n_states; ++i)
            {
                norm2_first_arg = norm2_first_arg +
                                  param_fn([&X, i, k]() { return -X(i, k); }) * var("X", {i, k});
            }

            // +(-u0^T)*u
            for (size_t i = 0; i < n_inputs; ++i)
            {
                norm2_first_arg = norm2_first_arg +
                                  param_fn([&U, i, k]() { return -U(i, k); }) * var("U", {i, k});
            }

            // +(-0.5)*Delta
            norm2_first_arg = norm2_first_arg + (-0.5) * var("Delta", {k});

            // +(0.5 + 0.5*x0^T*x0 + 0.5*u0^T*u0)
            norm2_first_arg = norm2_first_arg + param_fn([&X, &U, k]() {
                                  return 0.5 * (1.0 + X.col(k).dot(X.col(k)) + U.col(k).dot(U.col(k)));
                              });

            norm2_args.push_back(norm2_first_arg);
        }

        // (I)*x,
        for (size_t i = 0; i < n_states; ++i)
        {
            norm2_args.push_back((1.0) * var("X", {i, k}));
        }

        // (I)*u
        for (size_t i = 0; i < n_inputs; ++i)
        {
            norm2_args.push_back((1.0) * var("U", {i, k}));
        }

        // Right hand side
        // ( x0^T)*x  +( u0^T)*u  +( 0.5)*Delta  +(0.5 - 0.5*x0^T*x0 - 0.5*u0^T*u0)
        optimization_problem::AffineExpression rhs;

        // ( x0^T)*x
        for (size_t i = 0; i < n_states; ++i)
        {
            rhs = rhs + param(X(i, k)) * var("X", {i, k});
        }

        // +( u0^T)*u
        for (size_t i = 0; i < n_inputs; ++i)
        {
            rhs = rhs + param(U(i, k)) * var("U", {i, k});
        }

        // +( 0.5)*Delta
        rhs = rhs + (0.5) * var("Delta", {k});

        // +(0.5 - 0.5*x0^T*x0 - 0.5*u0^T*u0)
        rhs = rhs + param_fn([&X, &U, k]() {
                  return 0.5 * (1.0 - X.col(k).dot(X.col(k)) - U.col(k).dot(U.col(k)));
              });

        socp.add_constraint(optimization_problem::norm2(norm2_args) <= rhs);
    }

    /*
     * Build combined state/input trust region over all K:
     *   norm2([ Delta(1), Delta(2), ... , Delta(K) ]) <= norm2_Delta
     */
    {
        vector<optimization_problem::AffineExpression> norm2_args;
        for (size_t k = 0; k < K; k++)
        {
            norm2_args.push_back((1.0) * var("Delta", {k}));
        }
        socp.add_constraint(optimization_problem::norm2(norm2_args) <= (1.0) * var("norm2_Delta", {}));

        // Minimize norm2_Delta
        socp.add_minimization_term(param(weight_trust_region_xu) * var("norm2_Delta", {}));
    }

    socp.add_constraint((1.0) * var("sigma", {}) + (-0.01) >= (0.0)); // Total time must not be negative
    model.add_application_constraints(socp, X, U);
    return socp;
}