#include "SCProblem.hpp"

namespace scpp
{

    std::shared_ptr<cvx::OptimizationProblem> buildSCProblem(
        double &weight_time,
        double &weight_trust_region_time,
        double &weight_trust_region_trajectory,
        double &weight_virtual_control,
        trajectory_data_t &td,
        discretization_data_t &dd)
    {
        const size_t K = td.n_X();

        auto socp = std::make_shared<cvx::OptimizationProblem>();

        cvx::MatrixX v_X = socp->addVariable("X", Model::state_dim, td.n_X());                   // states
        cvx::MatrixX v_U = socp->addVariable("U", Model::input_dim, td.n_U());                   // inputs
        cvx::MatrixX v_nu = socp->addVariable("nu", Model::state_dim, td.n_X() - 1);             // virtual control
        cvx::MatrixX v_nu_bound = socp->addVariable("nu_bound", Model::state_dim, td.n_X() - 1); // virtual control
        cvx::Scalar v_norm1_nu = socp->addVariable("norm1_nu");                                  // virtual control norm upper bound
        cvx::VectorX v_delta = socp->addVariable("delta", K);                                    // change of the stacked [ x(k), u(k) ] vector

        cvx::Scalar v_sigma;
        cvx::Scalar v_delta_sigma;
        if (dd.variableTime())
        {
            v_sigma = socp->addVariable("sigma");
            v_delta_sigma = socp->addVariable("delta_sigma"); // squared change of sigma
            // minimize total time
            socp->addCostTerm(cvx::dynpar(weight_time) * v_sigma);
            // Total time must not be negative
            socp->addConstraint(cvx::greaterThan(v_sigma, 0.001));
        }

        for (size_t k = 0; k < K - 1; k++)
        {
            /**
             * Build linearized model equality constraint
             *    x(k+1) == A x(k) + B u(k) + C u(k+1) + Sigma sigma + z + nu
             * 
             */
            cvx::VectorX lhs = cvx::dynpar(dd.A.at(k)) * v_X.col(k) +
                               cvx::dynpar(dd.B.at(k)) * v_U.col(k) +
                               cvx::dynpar(dd.z.at(k)) +
                               v_nu.col(k);

            if (dd.interpolatedInput())
            {
                lhs += cvx::dynpar(dd.C.at(k)) * v_U.col(k + 1);
            }
            if (dd.variableTime())
            {
                lhs += cvx::dynpar(dd.s.at(k)) * v_sigma;
            }

            socp->addConstraint(cvx::equalTo(lhs, v_X.col(k + 1)));
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
            socp->addConstraint(cvx::box(-v_nu_bound, v_nu, v_nu_bound));

            // sum(nu_bound) <= norm1_nu
            socp->addConstraint(cvx::lessThan(v_nu_bound.sum(), v_norm1_nu));

            // Minimize the virtual control
            socp->addCostTerm(cvx::dynpar(weight_virtual_control) * v_norm1_nu);
        }

        if (dd.variableTime())
        {
            /**
             *  Build sigma trust region
             * (sigma - sigma0) * (sigma - sigma0) <= delta_sigma
             *          is equivalent to
             * norm2(
             *        0.5 - 0.5 * delta_sigma
             *        sigma0 - sigma
             *      )
             *      <= 0.5 + 0.5 * delta_sigma;
             */
            {
                cvx::VectorX norm2_terms(2);
                norm2_terms << cvx::par(0.5) + cvx::par(-0.5) * v_delta_sigma,
                    -cvx::dynpar(td.t) + v_sigma;

                socp->addConstraint(cvx::lessThan(norm2_terms.norm(), cvx::par(0.5) + cvx::par(0.5) * v_delta_sigma));

                // Minimize delta_sigma
                socp->addCostTerm(cvx::dynpar(weight_trust_region_time) * v_delta_sigma);
            }
        }

        for (size_t k = 0; k < K; k++)
        {
            /**
             * Build state and input trust-region:
             *
             *  norm2(
             *        (x - x0)
             *        (u - u0)
             *      )
             *     <= delta;
             *
             */

            cvx::VectorX norm2_terms(v_X.rows());
            norm2_terms << cvx::dynpar(td.X.at(k)) - v_X.col(k);

            if (dd.interpolatedInput() or (not dd.interpolatedInput() and k < K - 1))
            {
                norm2_terms.conservativeResize(v_X.rows() + v_U.rows());
                norm2_terms.tail(v_U.rows()) = cvx::dynpar(td.U.at(k)) - v_U.col(k);
            }

            socp->addConstraint(cvx::lessThan(norm2_terms.norm(), v_delta(k)));
        }

        /**
         * Minimize combined state/input trust region over all K:
         * 
         */
        {
            // Minimize trust region cost
            socp->addCostTerm(cvx::dynpar(weight_trust_region_trajectory) * v_delta.sum());
        }

        return socp;
    }

} // namespace scpp
