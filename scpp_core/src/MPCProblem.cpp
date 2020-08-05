#include "MPCProblem.hpp"

namespace scpp
{

cvx::OptimizationProblem buildMPCProblem(
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
    cvx::OptimizationProblem socp;

    cvx::MatrixX v_X = socp.addVariable("X", Model::state_dim, X.size()); // states
    cvx::MatrixX v_U = socp.addVariable("U", Model::input_dim, U.size()); // inputs
    cvx::Scalar v_error_cost = socp.addVariable("error_cost");           // error minimization term
    cvx::Scalar v_input_cost = socp.addVariable("input_cost");           // input minimization term

    // Initial state
    for (size_t i = 0; i < Model::state_dim; i++)
    {
        socp.addConstraint(cvx::equalTo(v_X.col(0), cvx::dynpar(x_init)));
    }

    for (size_t k = 0; k < X.size() - 1; k++)
    {
        /**
         * Build linearized model equality constraint
         *    x(k+1) == A x(k) + B u(k) + z
         * 
         */
        cvx::VectorX lhs;
        if (constant_dynamics)
        {
            lhs = cvx::par(A) * v_X.col(k) +
                  cvx::par(B) * v_U.col(k) +
                  cvx::par(z);
        }
        else
        {
            lhs = cvx::dynpar(A) * v_X.col(k) +
                  cvx::dynpar(B) * v_U.col(k) +
                  cvx::dynpar(z);
        }

        socp.addConstraint(cvx::equalTo(lhs, v_X.col(k + 1)));
    }

    /**
     * Build error cost
     * 
     */
    cvx::VectorX error_norm2_args(1 * v_X.rows());
    if (intermediate_cost_active)
    {
        error_norm2_args.resize((X.size()-1) * v_X.rows());
        for (size_t k = 1; k < X.size() - 1; k++)
        {
            error_norm2_args.segment((k-1) * v_X.cols(), v_X.cols()) = cvx::dynpar(state_weights_intermediate).cwiseProduct(-cvx::dynpar(x_final) + v_X.col(k));
        }
    }
    error_norm2_args.tail(v_X.rows()) = cvx::dynpar(state_weights_terminal).cwiseProduct(-cvx::dynpar(x_final) + v_X.col(v_X.cols() - 1));
    socp.addConstraint(cvx::lessThan(error_norm2_args.norm(), v_error_cost));
    socp.addCostTerm(v_error_cost);

    /**
     * Build input cost
     * 
     */
    cvx::VectorX input_norm2_args(U.size() * v_U.rows());
    for (size_t k = 0; k < U.size(); k++)
    {
        input_norm2_args.segment(k * v_U.rows(), v_U.rows()) = cvx::dynpar(input_weights).cwiseProduct(v_U.col(k));
    }
    socp.addConstraint(cvx::lessThan(input_norm2_args.norm(), v_input_cost));
    socp.addCostTerm(v_input_cost);

    return socp;
}

} // namespace scpp
