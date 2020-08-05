#include "SCvxProblem.hpp"

namespace scpp
{

cvx::OptimizationProblem buildSCvxProblem(
    double &trust_region,
    double &weight_virtual_control,
    trajectory_data_t &td,
    discretization_data_t &dd)
{
    cvx::OptimizationProblem socp;

    cvx::MatrixX v_X = socp.addVariable("X", Model::state_dim, td.n_X());                   // states
    cvx::MatrixX v_U = socp.addVariable("U", Model::input_dim, td.n_U());                   // inputs
    cvx::MatrixX v_nu = socp.addVariable("nu", Model::state_dim, td.n_X() - 1);             // virtual control
    cvx::MatrixX v_nu_bound = socp.addVariable("nu_bound", Model::state_dim, td.n_X() - 1); // virtual control lower/upper bound
    cvx::Scalar v_norm1_nu = socp.addVariable("norm1_nu");                                 // virtual control norm

    for (size_t k = 0; k < td.n_X() - 1; k++)
    {
        /**
         * Build linearized model equality constraint
         *    x(k+1) == A x(k) + B u(k) + C u(k+1) + z + nu
         * 
         */
        cvx::VectorX lhs = cvx::dynpar(dd.A.at(k)) * v_X.col(k) +
                           cvx::dynpar(dd.B.at(k)) * v_U.col(k) +
                           cvx::dynpar(dd.z.at(k)) +
                           v_nu.col(k);

        if (td.interpolatedInput())
        {
            lhs += cvx::dynpar(dd.C.at(k)) * v_U.col(k + 1);
        }

        socp.addConstraint(cvx::equalTo(lhs, v_X.col(k + 1)));
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
        socp.addConstraint(cvx::box(-v_nu_bound, v_nu, v_nu_bound));

        // sum(nu_bound) <= norm1_nu
        socp.addConstraint(cvx::lessThan(v_nu_bound.sum(), v_norm1_nu));

        // Minimize the virtual control
        socp.addCostTerm(cvx::dynpar(weight_virtual_control) * v_norm1_nu);
    }

    for (size_t k = 0; k < td.n_U(); k++)
    {
        /**
         * Build input trust region:
         *     norm2(u - u0)  <=  trust_region
         *
         */
        cvx::VectorX norm2_args = cvx::dynpar(td.U.at(k)) + -v_U.col(k);

        socp.addConstraint(cvx::lessThan(norm2_args.norm(), cvx::dynpar(trust_region)));
    }

    return socp;
}

} // namespace scpp
