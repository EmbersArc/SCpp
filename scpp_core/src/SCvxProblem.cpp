#include "optimizationProblem.hpp"
#include "SCvxProblem.hpp"

namespace scpp
{

op::SecondOrderConeProgram buildSCvxProblem(
    double &trust_region,
    double &weight_virtual_control,
    trajectory_data_t &td,
    discretization_data_t &dd)
{
    op::SecondOrderConeProgram socp;

    op::Variable v_X = socp.createVariable("X", Model::state_dim, td.n_X());                   // states
    op::Variable v_U = socp.createVariable("U", Model::input_dim, td.n_U());                   // inputs
    op::Variable v_nu = socp.createVariable("nu", Model::state_dim, td.n_X() - 1);             // virtual control
    op::Variable v_nu_bound = socp.createVariable("nu_bound", Model::state_dim, td.n_X() - 1); // virtual control lower/upper bound
    op::Variable v_norm1_nu = socp.createVariable("norm1_nu");                                 // virtual control norm

    for (size_t k = 0; k < td.n_X() - 1; k++)
    {
        /**
         * Build linearized model equality constraint
         *    x(k+1) == A x(k) + B u(k) + C u(k+1) + z + nu
         *   -x(k+1)  + A x(k) + B u(k) + C u(k+1) + z + nu == 0
         * 
         */
        socp.addConstraint(op::Parameter(-1.0) * v_X.col(k + 1) +
                               op::Parameter(&dd.A.at(k)) * v_X.col(k) +
                               op::Parameter(&dd.B.at(k)) * v_U.col(k) +
                               op::Parameter(&dd.C.at(k)) * v_U.col(k + 1) +
                               op::Parameter(&dd.z.at(k)) +
                               op::Parameter(1.0) * v_nu.col(k) ==
                           0.);
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
        socp.addConstraint(op::Parameter(1.0) * v_nu_bound + op::Parameter(1.0) * v_nu >= 0.);
        socp.addConstraint(op::Parameter(1.0) * v_nu_bound + op::Parameter(-1.0) * v_nu >= 0.);

        op::AffineExpression bound_sum;
        for (size_t row = 0; row < v_nu.rows(); row++)
        {
            for (size_t col = 0; col < v_nu.cols(); col++)
            {
                bound_sum = bound_sum + op::ParameterSource(-1.0) * v_nu_bound(row, col);
            }
        }

        // sum(nu_bound) <= norm1_nu
        socp.addConstraint(op::Parameter(1.0) * v_norm1_nu + op::Affine(bound_sum) >= (0.0));

        // Minimize the virtual control
        socp.addMinimizationTerm(op::Parameter(&weight_virtual_control) * v_norm1_nu);
    }

    for (size_t k = 0; k < td.n_U(); k++)
    {
        /**
         * Build input trust region:
         *     norm2(u - u0)  <=  trust_region
         *
         */
        op::Affine norm2_args = op::Parameter(&td.U.at(k)) + op::Parameter(-1.0) * v_U.col(k);

        socp.addConstraint(op::Norm2(norm2_args) <= op::Parameter(&trust_region));
    }

    return socp;
}

} // namespace scpp
