#include "rocket2d.hpp"
#include "common.hpp"

namespace scpp::models
{

    void Rocket2d::systemFlowMap(const state_vector_ad_t &x,
                                 const input_vector_ad_t &u,
                                 const param_vector_ad_t &par,
                                 state_vector_ad_t &f)
    {
        using T = scalar_ad_t;

        auto m = par(0);
        auto J_B = par(1);
        auto g_I = par.segment<2>(2);
        auto r_T_B = par.segment<2>(4);
        // = 6 parameters

        // state variables
        // x, y, vx, vy, eta, omega
        Eigen::Matrix<T, 2, 1> v = x.segment<2>(2);
        T eta = x(4);
        T w = x(5);

        // input variables
        T angle = u(0);
        T magnitude = u(1);
        Eigen::Matrix<T, 2, 1> T_B =
            Eigen::Rotation2D<T>(angle) * Eigen::Matrix<T, 2, 1>(0., magnitude);

        Eigen::Rotation2D<T> R_I_B(eta);

        f.segment<2>(0) << v;
        f.segment<2>(2) << 1. / m * (R_I_B * T_B) + g_I;
        f(4) = w;
        f(5) = 1. / J_B * (r_T_B.x() * T_B.y() - r_T_B.y() * T_B.x());
    }

    void Rocket2d::getOperatingPoint(state_vector_t &x, input_vector_t &u)
    {
        x.setZero();
        u << 0, -p.g_I * p.m;
    }

    void Rocket2d::addApplicationConstraints(std::shared_ptr<cvx::OptimizationProblem> socp,
                                             state_vector_v_t &,
                                             input_vector_v_t &)
    {
        cvx::MatrixX v_X, v_U;
        socp->getVariable("X", v_X);
        socp->getVariable("U", v_U);

        if (p.constrain_initial_final)
        {
            // Initial and final state
            socp->addConstraint(cvx::equalTo(cvx::dynpar(p.x_init), v_X.col(0)));
            socp->addConstraint(cvx::equalTo(cvx::dynpar(p.x_final), v_X.rightCols(1)));
            socp->addConstraint(cvx::equalTo(v_U(0, v_U.cols() - 1), 0.));
        }

        // State Constraints:
        // Glideslope
        socp->addConstraint(cvx::lessThan(v_X.row(0).colwise().norm(),
                                          cvx::dynpar(p.tan_gamma_gs) * v_X.row(1)));
        // // Max Tilt Angle
        socp->addConstraint(cvx::box(-cvx::dynpar(p.theta_max),
                                     v_X.row(4),
                                     cvx::dynpar(p.theta_max)));
        // Max Rotation Velocity
        socp->addConstraint(cvx::box(-cvx::dynpar(p.w_B_max),
                                     v_X.row(5),
                                     cvx::dynpar(p.w_B_max)));

        // Control Constraints:
        // Gimbal Range
        socp->addConstraint(cvx::box(-cvx::dynpar(p.gimbal_max),
                                     v_U.row(0),
                                     cvx::dynpar(p.gimbal_max)));
        // Thrust Range
        socp->addConstraint(cvx::box(cvx::dynpar(p.T_min),
                                     v_U.row(1),
                                     cvx::dynpar(p.T_max)));
    }

    void Rocket2d::nondimensionalize()
    {
        p.nondimensionalize();
    }

    void Rocket2d::redimensionalize()
    {
        p.redimensionalize();
    }

    void Rocket2d::nondimensionalizeTrajectory(trajectory_data_t &td)
    {
        for (auto &x : td.X)
        {
            x.head(4) /= p.r_scale;
        }
        for (auto &u : td.U)
        {
            u(1) /= p.m_scale * p.r_scale;
        }
    }

    void Rocket2d::redimensionalizeTrajectory(trajectory_data_t &td)
    {
        for (auto &x : td.X)
        {
            x.head(4) *= p.r_scale;
        }
        for (auto &u : td.U)
        {
            u(1) *= p.m_scale * p.r_scale;
        }
    }

    void Rocket2d::getInitializedTrajectory(trajectory_data_t &td)
    {
        for (size_t k = 0; k < td.n_X(); k++)
        {
            const double alpha1 = double(td.n_X() - k) / td.n_X();
            const double alpha2 = double(k) / td.n_X();

            td.X.at(k) = alpha1 * p.x_init + alpha2 * p.x_final;
        }

        for (auto &u : td.U)
        {
            u << 0, (p.T_max + p.T_min) / 2;
        }

        td.t = p.final_time;
    }

    void Rocket2d::loadParameters()
    {
        p.loadFromFile(getParameterFolder() + "/model.info");
    }

    void Rocket2d::getNewModelParameters(param_vector_t &par)
    {
        par << p.m, p.J_B, p.g_I, p.r_T_B;

        p.tan_gamma_gs = std::tan(p.gamma_gs);
    }

    void Rocket2d::Parameters::loadFromFile(const std::string &path)
    {
        ParameterServer param(path);

        param.loadMatrix("g_I", g_I);
        param.loadScalar("J_B", J_B);
        param.loadMatrix("r_T_B", r_T_B);

        Eigen::Vector2d r_init, v_init;
        Eigen::Vector2d r_final, v_final;
        double w_init, w_final;

        param.loadMatrix("r_init", r_init);
        param.loadMatrix("v_init", v_init);
        param.loadScalar("eta_init", eta_init);
        param.loadScalar("w_init", w_init);

        param.loadMatrix("r_final", r_final);
        param.loadMatrix("v_final", v_final);
        param.loadScalar("eta_final", eta_final);
        param.loadScalar("w_final", w_final);

        param.loadScalar("final_time", final_time);

        param.loadScalar("m", m);
        param.loadScalar("T_min", T_min);
        param.loadScalar("T_max", T_max);
        param.loadScalar("gamma_gs", gamma_gs);
        param.loadScalar("gimbal_max", gimbal_max);
        param.loadScalar("theta_max", theta_max);
        param.loadScalar("w_B_max", w_B_max);
        param.loadScalar("constrain_initial_final", constrain_initial_final);
        param.loadScalar("add_slack_variables", add_slack_variables);

        deg2rad(gimbal_max);
        deg2rad(theta_max);
        deg2rad(gamma_gs);
        deg2rad(w_B_max);
        deg2rad(w_init);
        deg2rad(w_final);
        deg2rad(eta_init);
        deg2rad(eta_final);

        x_init << r_init, v_init, eta_init, w_init;
        x_final << r_final, v_final, eta_final, w_final;
    }

    void Rocket2d::Parameters::nondimensionalize()
    {
        r_scale = x_init.head(2).norm();
        m_scale = m;

        m /= m_scale;
        r_T_B /= r_scale;
        g_I /= r_scale;
        J_B /= m_scale * r_scale * r_scale;

        x_init.segment<2>(0) /= r_scale;
        x_init.segment<2>(2) /= r_scale;

        x_final.segment<2>(0) /= r_scale;
        x_final.segment<2>(2) /= r_scale;

        T_min /= m_scale * r_scale;
        T_max /= m_scale * r_scale;
    }

    void Rocket2d::Parameters::redimensionalize()
    {
        m *= m_scale;
        r_T_B *= r_scale;
        g_I *= r_scale;
        J_B *= m_scale * r_scale * r_scale;

        x_init.segment<2>(0) *= r_scale;
        x_init.segment<2>(2) *= r_scale;

        x_final.segment<2>(0) *= r_scale;
        x_final.segment<2>(2) *= r_scale;

        T_min *= m_scale * r_scale;
        T_max *= m_scale * r_scale;
    }

} // namespace scpp::models