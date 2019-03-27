#include "common.hpp"
#include "cartpole.hpp"

using std::string;
using std::vector;

namespace cartpole
{

Cartpole::Cartpole()
{
    p.loadFromFile();
}

void Cartpole::systemFlowMap(const state_vector_ad_t &x,
                             const input_vector_ad_t &u,
                             const param_vector_ad_t &p,
                             state_vector_ad_t &f)
{
    typedef scalar_ad_t T;

    // X = [x, dx, theta, dtheta]^T
    T g = T(this->p.g);
    T l = T(this->p.l);
    T m_c = T(this->p.m_c);
    T m_p = T(this->p.m_p);

    T num = -g * sin(x(2)) + cos(x(2)) * (-u(0) - m_p * l * x(3) * x(3) * sin(x(2))) / (m_c + m_p);
    T denum = l * (T(4. / 3.) - (m_p * cos(x(2)) * cos(x(2))) / (m_c + m_p));
    T ddtheta = num / denum;

    f(0) = x(1);
    f(1) = (u(0) + m_p * l * (x(3) * x(3) * sin(x(2)) - ddtheta * cos(x(2)))) / (m_c + m_p);
    f(2) = x(3);
    f(3) = ddtheta;
}

void Cartpole::getInitializedTrajectory(Eigen::MatrixXd &X,
                                        Eigen::MatrixXd &U)
{
    const size_t K = X.cols();

    for (size_t k = 0; k < K; k++)
    {
        const double alpha1 = double(K - k) / K;
        const double alpha2 = double(k) / K;

        X.col(k) = alpha1 * p.x_init + alpha2 * p.x_final;
        U(k) = sin(2. * M_PI * alpha2) * p.F_max / 4;
    }
}

void Cartpole::addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                         Eigen::MatrixXd &X0,
                                         Eigen::MatrixXd &U0)
{
    const size_t K = X0.cols();

    auto var = [&socp](const string &name, const vector<size_t> &indices = {}) { return socp.getVariable(name, indices); };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    // auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    // Initial state
    socp.addConstraint((-1.0) * var("X", {0, 0}) + param(p.x_init(0)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {1, 0}) + param(p.x_init(1)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {2, 0}) + param(p.x_init(2)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {3, 0}) + param(p.x_init(3)) == 0.0);

    // Final state
    socp.addConstraint((-1.0) * var("X", {0, K - 1}) + param(p.x_final(0)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {1, K - 1}) + param(p.x_final(1)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {2, K - 1}) + param(p.x_final(2)) == 0.0);
    socp.addConstraint((-1.0) * var("X", {3, K - 1}) + param(p.x_final(3)) == 0.0);

    // Input
    for (size_t k = 0; k < K; k++)
    {
        socp.addConstraint((1.0) * var("U", {0, k}) + (p.F_max) >= (0.0));
        socp.addConstraint((-1.0) * var("U", {0, k}) + (p.F_max) >= (0.0));
    }
    // State
    for (size_t k = 0; k < K; k++)
    {
        socp.addConstraint((1.0) * var("X", {0, k}) + (p.x_max) >= (0.0));
        socp.addConstraint((-1.0) * var("X", {0, k}) + (p.x_max) >= (0.0));

        socp.addConstraint((1.0) * var("X", {3, k}) + (p.dtheta_max) >= (0.0));
        socp.addConstraint((-1.0) * var("X", {3, k}) + (p.dtheta_max) >= (0.0));
    }
}

void Cartpole::nondimensionalize()
{
    p.nondimensionalize();
}

void Cartpole::redimensionalize()
{
    p.redimensionalize();
}

void Cartpole::nondimensionalizeTrajectory(Eigen::MatrixXd &X,
                                           Eigen::MatrixXd &U)
{
    p.nondimensionalizeTrajectory(X, U);
}

void Cartpole::redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                          Eigen::MatrixXd &U)
{
    p.redimensionalizeTrajectory(X, U);
}

void Cartpole::getNewModelParameters(param_vector_t &param)
{
}

void Cartpole::Parameters::loadFromFile()
{
    ParameterServer param(fmt::format("../models/config/{}.info", getModelName()));

    param.loadMatrix("x_init", x_init);
    param.loadMatrix("x_final", x_final);
    param.loadScalar("m_c", m_c);
    param.loadScalar("m_p", m_p);
    param.loadScalar("l", l);
    param.loadScalar("g", g);

    param.loadScalar("x_max", x_max);
    param.loadScalar("F_max", F_max);
    param.loadScalar("dtheta_max", dtheta_max);

    deg2rad(dtheta_max);
    deg2rad(x_init(2));
    deg2rad(x_init(3));
    deg2rad(x_final(2));
    deg2rad(x_final(3));
}

void Cartpole::Parameters::nondimensionalize()
{
    m_scale = m_c + m_p;
    r_scale = l;

    x_init.head(2) /= r_scale;
    x_final.head(2) /= r_scale;

    m_c /= m_scale;
    m_p /= m_scale;
    l /= r_scale;
    g /= r_scale;

    x_max /= r_scale;
    F_max /= m_scale * r_scale;
}

void Cartpole::Parameters::redimensionalize()
{
    x_init.head(2) *= r_scale;
    x_final.head(2) *= r_scale;

    m_c *= m_scale;
    m_p *= m_scale;
    l *= r_scale;
    g *= r_scale;

    x_max *= r_scale;
    F_max *= m_scale * r_scale;
}

void Cartpole::Parameters::nondimensionalizeTrajectory(Eigen::MatrixXd &X, Eigen::MatrixXd &U) const
{
    X.topRows(2) /= r_scale;
    U /= m_scale * r_scale;
}

void Cartpole::Parameters::redimensionalizeTrajectory(Eigen::MatrixXd &X, Eigen::MatrixXd &U) const
{
    X.topRows(2) *= r_scale;
    U *= m_scale * r_scale;
}

} // namespace cartpole