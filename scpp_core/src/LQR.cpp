#include "LQR.hpp"

constexpr size_t STATE_DIM = Model::state_dim;
using schur_matrix_t = Eigen::Matrix<double, 2 * STATE_DIM, 2 * STATE_DIM>;
using factor_matrix_t = Eigen::Matrix<double, 2 * STATE_DIM, STATE_DIM>;

bool solveSchurIterative(const schur_matrix_t &M,
                         Model::state_matrix_t &P,
                         double epsilon,
                         size_t maxIterations)
{
    bool converged = false;

    schur_matrix_t Mlocal = M;

    size_t iterations = 0;
    while (not converged)
    {
        if (iterations > maxIterations)
            return false;

        const schur_matrix_t Mdiff = Mlocal - Mlocal.inverse();

        const schur_matrix_t Mnew = Mlocal - 0.5 * Mdiff;

        converged = Mnew.isApprox(Mlocal, epsilon);

        Mlocal = Mnew;

        iterations++;
    }

    /* break down M and extract M11 M12 M21 M22 */
    const Model::state_matrix_t M11(Mlocal.topLeftCorner<STATE_DIM, STATE_DIM>());
    const Model::state_matrix_t M12(Mlocal.topRightCorner<STATE_DIM, STATE_DIM>());
    const Model::state_matrix_t M21(Mlocal.bottomLeftCorner<STATE_DIM, STATE_DIM>());
    const Model::state_matrix_t M22(Mlocal.bottomRightCorner<STATE_DIM, STATE_DIM>());

    /* find M and N using the elements of M	*/
    factor_matrix_t U;
    factor_matrix_t V;

    U.topRows<STATE_DIM>() = M12;
    U.bottomRows<STATE_DIM>() = M22 + Model::state_matrix_t::Identity();

    V.topRows<STATE_DIM>() = M11 + Model::state_matrix_t::Identity();
    V.bottomRows<STATE_DIM>() = M21;

    /* Solve for S from the equation   MS=N */
    Eigen::FullPivLU<factor_matrix_t> FullPivLU_;
    FullPivLU_.compute(U);

    P = FullPivLU_.solve(-V);

    return true;
}

bool careSolve(const Model::state_matrix_t &Q,
               const Model::input_matrix_t &R,
               const Model::state_matrix_t &A,
               const Model::control_matrix_t &B,
               Model::state_matrix_t &P,
               Model::input_matrix_t &R_inverse)
{
    if ((R - Model::input_matrix_t::Identity().cwiseProduct(R)).any())
    {
        R_inverse = R.inverse();
    }
    else
    {
        R_inverse.setZero();
        R_inverse.diagonal() = R.diagonal().cwiseInverse();
    }

    schur_matrix_t M;
    M << A, -B * R_inverse * B.transpose(), -Q, -A.transpose();

    return solveSchurIterative(M, P, 1e-8, 100);
}

bool ComputeLQR(const Model::state_matrix_t &Q,
                const Model::input_matrix_t &R,
                const Model::state_matrix_t &A,
                const Model::control_matrix_t &B,
                Model::feedback_matrix_t &K)
{
    // fmt::print("[LQR] Checking controllability...\n");

    {
        using ctrl_matrix_t = Eigen::Matrix<double, Model::state_dim, Model::state_dim * Model::input_dim>;
        ctrl_matrix_t C;
        C.block<Model::state_dim, Model::input_dim>(0, 0) = B;
        for (size_t i = 1; i < Model::state_dim; i++)
        {
            C.block<Model::state_dim, Model::input_dim>(0, i * Model::input_dim).noalias() =
                A * C.block<Model::state_dim, Model::input_dim>(0, (i - 1) * Model::input_dim);
        }

        // check if system is controllable
        assert(Eigen::FullPivLU<ctrl_matrix_t>(C).rank() == Model::state_dim);
    }

    // fmt::print("[LQR] Computing feedback gains.\n");
    Model::input_matrix_t R_inverse;
    Model::state_matrix_t P;
    bool success = careSolve(Q, R, A, B, P, R_inverse);
    K = (R_inverse * (B.transpose() * P));

    assert(not K.hasNaN());

    return success;
}
