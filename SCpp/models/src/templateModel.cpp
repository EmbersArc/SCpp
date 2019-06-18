#include "templateModel.hpp"

using std::string;
using std::vector;

namespace templateModel
{

TemplateModel::TemplateModel()
{
}

void TemplateModel::systemFlowMap(const state_vector_ad_t &x,
                                  const input_vector_ad_t &u,
                                  const param_vector_ad_t &p,
                                  state_vector_ad_t &f)
{
    typedef scalar_ad_t T;

    f(0) = x(0) + u(0);
}

void TemplateModel::getInitializedTrajectory(Eigen::MatrixXd &X,
                                             Eigen::MatrixXd &U)
{
    const size_t K = X.cols();

    for (size_t k = 0; k < K; k++)
    {
        const double alpha1 = double(K - k) / K;
        const double alpha2 = double(k) / K;

        X.col(k) = alpha1 * p.x_init + alpha2 * p.x_final;
    }
}

void TemplateModel::addApplicationConstraints(op::SecondOrderConeProgram &socp,
                                              Eigen::MatrixXd &X0,
                                              Eigen::MatrixXd &U0)
{
    const size_t K = X0.cols();

    auto var = [&socp](const string &name, const vector<size_t> &indices = {}) { return socp.getVariable(name, indices); };
    auto param = [](double &param_value) { return op::Parameter(&param_value); };
    auto param_fn = [](std::function<double()> callback) { return op::Parameter(callback); };

    // Initial state
    for (size_t i = 0; i < STATE_DIM_; i++)
    {
        socp.addConstraint((-1.0) * var("X", {i, 0}) + param(p.x_init(i)) == 0.0);
    }
    // Final state
    for (size_t i = 0; i < STATE_DIM_; i++)
    {
        socp.addConstraint((-1.0) * var("X", {i, 0}) + param(p.x_final(i)) == 0.0);
    }
}

void TemplateModel::nondimensionalize()
{
    p.nondimensionalize();
}

void TemplateModel::redimensionalize()
{
    p.redimensionalize();
}

void TemplateModel::nondimensionalizeTrajectory(Eigen::MatrixXd &X,
                                                Eigen::MatrixXd &U)
{
    p.nondimensionalizeTrajectory(X, U);
}

void TemplateModel::redimensionalizeTrajectory(Eigen::MatrixXd &X,
                                               Eigen::MatrixXd &U)
{
    p.redimensionalizeTrajectory(X, U);
}

void TemplateModel::getNewModelParameters(param_vector_t &param)
{
    param << p.m;
}

void TemplateModel::Parameters::loadFromFile()
{
    ParameterServer param(fmt::format("../SCpp/models/config/{}.info", getModelName()));

    Eigen::Vector3d r_init;
    Eigen::Vector3d r_final;

    param.loadMatrix("r_init", r_init);
    param.loadMatrix("r_final", r_final);

    x_init << r_init;
    x_final << r_final;
}

void TemplateModel::Parameters::nondimensionalize()
{
    m_scale = 1;
    r_scale = 1;
}

void TemplateModel::Parameters::redimensionalize()
{
}

void TemplateModel::Parameters::nondimensionalizeTrajectory(Eigen::MatrixXd &X, Eigen::MatrixXd &U) const
{
}

void TemplateModel::Parameters::redimensionalizeTrajectory(Eigen::MatrixXd &X, Eigen::MatrixXd &U) const
{
}

} // namespace templateModel