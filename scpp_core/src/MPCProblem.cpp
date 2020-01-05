#include "optimizationProblem.hpp"
#include "MPCProblem.hpp"

namespace scpp
{

op::SecondOrderConeProgram buildMPCProblem(
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
    op::SecondOrderConeProgram socp;

    op::Variable v_X = socp.createVariable("X", Model::state_dim, X.size()); // states
    op::Variable v_U = socp.createVariable("U", Model::input_dim, U.size()); // inputs
    op::Variable v_error_cost = socp.createVariable("error_cost");           // error minimization term
    op::Variable v_input_cost = socp.createVariable("input_cost");           // input minimization term

    // Initial state
    for (size_t i = 0; i < Model::state_dim; i++)
    {
        socp.addConstraint(v_X.col(0) == op::Parameter(&x_init));
    }

    for (size_t k = 0; k < X.size() - 1; k++)
    {
        /**
         * Build linearized model equality constraint
         *    x(k+1) == A x(k) + B u(k) + z
         * 
         */
        op::Affine lhs;
        if (constant_dynamics)
        {
            lhs = op::Parameter(A) * v_X.col(k) +
                  op::Parameter(B) * v_U.col(k) +
                  op::Parameter(z);
        }
        else
        {
            lhs = op::Parameter(&A) * v_X.col(k) +
                  op::Parameter(&B) * v_U.col(k) +
                  op::Parameter(&z);
        }

        socp.addConstraint(lhs == v_X.col(k + 1));
    }

    /**
     * Build error cost
     * 
     */
    op::Affine error_norm2_args;
    if (intermediate_cost_active)
    {
        for (size_t k = 1; k < X.size() - 1; k++)
        {
            error_norm2_args = op::vstack({error_norm2_args,
                                           op::Parameter(&state_weights_intermediate).cwiseProduct(-op::Parameter(&x_final) + v_X.col(k))});
        }
    }
    error_norm2_args = op::vstack({error_norm2_args,
                                   op::Parameter(&state_weights_terminal).cwiseProduct(-op::Parameter(&x_final) + v_X.col(v_X.cols() - 1))});
    socp.addConstraint(op::Norm2(error_norm2_args) <= v_error_cost);
    socp.addMinimizationTerm(v_error_cost);

    /**
     * Build input cost
     * 
     */
    op::Affine input_norm2_args;
    for (size_t k = 0; k < U.size(); k++)
    {
        input_norm2_args = op::vstack({input_norm2_args,
                                       op::Parameter(&input_weights).cwiseProduct(v_U.col(k))});
    }
    socp.addConstraint(op::Norm2(input_norm2_args) <= v_input_cost);
    socp.addMinimizationTerm(v_input_cost);

    return socp;
}

} // namespace scpp
