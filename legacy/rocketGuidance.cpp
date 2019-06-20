#include "rocketGuidance.hpp"

using std::get;
namespace chrono = std::chrono;

std::mutex g_mutex;

RocketGuidance::RocketGuidance()
    : m_connection(krpc::connect()), m_model(std::make_shared<Model>())
{
    m_krpc = std::make_shared<KRPC>(&m_connection);
    m_spaceCenter = std::make_shared<SpaceCenter>(&m_connection);
    m_drawing = std::make_shared<Drawing>(&m_connection);

    m_solver = std::make_shared<FreeFinalTimeAlgorithm>(m_model);
    m_solver->initialize();

    m_vessel = m_spaceCenter->active_vessel();
    m_body = m_vessel.orbit().body();

    initializeReferenceFrames();

    m_warm_start = false;

    m_m_stream = m_vessel.mass_stream();
    m_r_I_stream = m_vessel.position_stream(m_ref_frame_pad);
    m_v_I_stream = m_vessel.velocity_stream(m_ref_frame_pad);
    m_q_I_stream = m_vessel.rotation_stream(m_ref_frame_pad);
    m_w_B_stream = m_vessel.angular_velocity_stream(m_ref_frame_nonrot);
    m_J_B_stream = m_vessel.moment_of_inertia_stream();
    m_engine = m_vessel.parts().engines()[0];
    m_r_T_B_stream = m_engine.thrusters()[0].thrust_position_stream(m_ref_frame_vessel);

    m_pitch_stream = m_vessel.control().pitch_stream();
    m_yaw_stream = m_vessel.control().yaw_stream();
    m_T_stream = m_vessel.control().throttle_stream();
}

void RocketGuidance::initializeReferenceFrames()
{
    // Define the landing site as launch pad
    double landing_latitude = -0.0972039;
    double landing_longitude = -74.5577;
    double landing_altitude = 18.1;

    // Determine landing site reference frame
    // (orientation: x=east, y=up, z=north)
    auto create_relative = [this](const SpaceCenter::ReferenceFrame &ref_frame, const tuple3d_t &location, const tuple4d_t &rotation = {0., 0., 0., 1.}) {
        return SpaceCenter::ReferenceFrame::create_relative(m_connection, ref_frame, location, rotation);
    };
    tuple3d_t landing_position = m_body.surface_position(landing_latitude, landing_longitude, m_body.reference_frame());
    Eigen::Quaterniond q_long;
    Eigen::Quaterniond q_lat;
    Eigen::Quaterniond q_rot;
    q_rot = Eigen::AngleAxisd(-M_PI_2, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(-M_PI_2, Eigen::Vector3d::UnitX());
    q_long = Eigen::AngleAxisd(-landing_longitude / 180. * M_PI, Eigen::Vector3d::UnitY());
    q_lat = Eigen::AngleAxisd(landing_latitude / 180. * M_PI, Eigen::Vector3d::UnitZ());
    Eigen::Quaterniond q_latlong = q_lat * q_long * q_rot;
    tuple4d_t q_latlong_tuple = std::make_tuple(q_latlong.x(), q_latlong.y(), q_latlong.z(), q_latlong.w());

    m_ref_frame_pad = create_relative(create_relative(m_body.reference_frame(), landing_position, q_latlong_tuple), std::make_tuple(0., landing_altitude, 0.));
    m_ref_frame_vessel = m_vessel.reference_frame();
    m_ref_frame_nonrot = m_body.non_rotating_reference_frame();
}

void RocketGuidance::drawTrajectory(const Eigen::MatrixXd &X)
{
    while (true)
    {
        m_drawing->clear();

        // coordinate frames
        auto dir_x_landing = m_drawing->add_line(std::make_tuple(0, 0, 0), std::make_tuple(30, 0, 0), m_ref_frame_pad, 10.);
        auto dir_y_landing = m_drawing->add_line(std::make_tuple(0, 0, 0), std::make_tuple(0, 30, 0), m_ref_frame_pad, 10.);
        auto dir_z_landing = m_drawing->add_line(std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 30), m_ref_frame_pad, 10.);
        dir_x_landing.set_color(std::make_tuple(1, 0, 0));
        dir_y_landing.set_color(std::make_tuple(0, 1, 0));
        dir_z_landing.set_color(std::make_tuple(0, 0, 1));
        auto dir_x_vessel = m_drawing->add_line(std::make_tuple(0, 0, 0), std::make_tuple(15, 0, 0), m_ref_frame_vessel, 10.);
        auto dir_y_vessel = m_drawing->add_line(std::make_tuple(0, 0, 0), std::make_tuple(0, 15, 0), m_ref_frame_vessel, 10.);
        auto dir_z_vessel = m_drawing->add_line(std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 15), m_ref_frame_vessel, 10.);
        dir_x_vessel.set_color(std::make_tuple(1, 0, 0));
        dir_y_vessel.set_color(std::make_tuple(0, 1, 0));
        dir_z_vessel.set_color(std::make_tuple(0, 0, 1));

        // trajectory and orientation
        size_t K = X.cols() - 1;
        for (size_t i = 0; i < K; i++)
        {
            tuple3d_t start = std::make_tuple(X(1, i), X(3, i), X(2, i));
            tuple3d_t end = std::make_tuple(X(1, i + 1), X(3, i + 1), X(2, i + 1));
            auto trajectory_line = m_drawing->add_line(start, end, m_ref_frame_pad);
            trajectory_line.set_thickness(1.);
            trajectory_line.set_color(std::make_tuple(1, 0.5, 0));

            Eigen::MatrixXd x = X.col(i);
            Eigen::Quaterniond q(x(7), x(8), x(9), x(10));
            Eigen::Vector3d attitude_rot = q.toRotationMatrix() * Eigen::Vector3d::UnitZ() * 10.;
            tuple3d_t end_attitude = std::make_tuple(get<0>(start) + attitude_rot(0), get<1>(start) + attitude_rot(2), get<2>(start) + attitude_rot(1));
            auto attitude_line = m_drawing->add_line(start, end_attitude, m_ref_frame_pad);
            attitude_line.set_thickness(1.);
        }

        std::this_thread::sleep_for(chrono::milliseconds(50));
    }
}

void RocketGuidance::prepareVessel()
{
    m_vessel.control().set_sas(false);
    m_vessel.control().set_rcs(false);
    m_vessel.control().set_reaction_wheels(false);
    m_engine.set_gimbal_locked(true);
    m_engine.set_gimbal_locked(false);
}

void RocketGuidance::updateModelParameters()
{
    m_model->p.final_time_guess = 10.;
    m_model->p.g_I << 0., 0., -m_body.surface_gravity();
    m_model->p.T_min = 0.2 * m_engine.max_thrust();
    m_model->p.T_max = m_engine.max_thrust();
    m_model->p.theta_max = 30. / 180. * M_PI;
    m_model->p.w_B_max = 20. / 180. * M_PI;
    const tuple3d_t r_T_B_current = m_r_T_B_stream();
    m_model->p.r_T_B << get<0>(r_T_B_current), get<2>(r_T_B_current), get<1>(r_T_B_current);
    const tuple3d_t m_J_B_current = m_J_B_stream();
    m_model->p.J_B << get<0>(m_J_B_current), get<2>(m_J_B_current), get<1>(m_J_B_current);
    m_model->p.alpha_m = 1. / (m_engine.specific_impulse() * std::abs(m_body.surface_gravity()));
    m_model->p.gimbal_max = m_engine.gimbal_range() * m_engine.gimbal_limit() / 180. * M_PI;
    m_model->p.gamma_gs = 45. / 180. * M_PI;

    const double m = m_m_stream();
    const tuple3d_t r_current = m_r_I_stream();
    const tuple3d_t v_current = m_v_I_stream();
    const tuple4d_t q_current = m_q_I_stream();
    const tuple3d_t w_current = m_spaceCenter->transform_direction(m_w_B_stream(), m_ref_frame_nonrot, m_ref_frame_vessel);
    const Eigen::Vector3d r(get<0>(r_current), get<2>(r_current), get<1>(r_current));
    const Eigen::Vector3d v(get<0>(v_current), get<2>(v_current), get<1>(v_current));
    const Eigen::Vector4d q(-get<3>(q_current), get<0>(q_current), get<2>(q_current), get<1>(q_current));
    const Eigen::Vector3d w(-get<0>(w_current), -get<2>(w_current), -get<1>(w_current));

    m_model->p.x_init << m, r, v, q, w;
    m_model->p.x_final << m_vessel.dry_mass(),
        0., 0., 0.,
        0., 0., 0.,
        1., 0., 0., 0.,
        0., 0., 0.;
}

void RocketGuidance::commandInput(const Eigen::MatrixXd &u)
{
    assert(u.rows() == 3 and u.cols() == 1);
    const double pitch = atan(u(1) / u(2)) / m_model->p.gimbal_max;
    const double yaw = -atan(u(0) / u(2)) / m_model->p.gimbal_max;
    m_vessel.control().set_pitch(pitch);
    m_vessel.control().set_yaw(yaw);
    m_vessel.control().set_throttle(u.norm() / m_model->p.T_max);
}

void RocketGuidance::getInterpolatedInput(Eigen::MatrixXd &u, const Eigen::MatrixXd &U, const double t)
{
    assert(0. <= t and t < 1.);
    const size_t K = U.cols() - 1;
    const double dt = 1. / K;
    const size_t k0 = size_t(t * K);
    const double delta = (t - k0 * dt) / dt;
    u = U.col(k0) + (U.col(k0 + 1) - U.col(k0)) * delta;
}

void RocketGuidance::recomputeTrajectory(Eigen::MatrixXd &X, Eigen::MatrixXd &U, double &t, bool &running)
{
    while (running)
    {
        updateModelParameters();
        chrono::time_point start_time = chrono::high_resolution_clock::now();
        m_solver->solve(m_warm_start.load());
        m_warm_start = true;
        {
            std::lock_guard<std::mutex> guard(g_mutex);
            m_start_time = start_time;
            m_solver->getSolution(X, U, t);
        }
        std::this_thread::sleep_for(chrono::milliseconds(1));
    }
}

void RocketGuidance::run()
{
    prepareVessel();

    Model::dynamic_matrix_t X, U;
    double t;

    bool running = true;
    m_solverThread = std::thread(&RocketGuidance::recomputeTrajectory, this, std::ref(X), std::ref(U), std::ref(t), std::ref(running));

    // wait for the first solution
    if (not m_warm_start)
    {
        m_krpc->set_paused(true);
        while (not m_warm_start)
        {
            std::this_thread::sleep_for(chrono::milliseconds(1));
            m_start_time = chrono::high_resolution_clock::now();
        }
        m_krpc->set_paused(false);
    }

    // running = false;
    // m_solverThread.join();
    m_drawingThread = std::thread(&RocketGuidance::drawTrajectory, this, std::ref(X));

    while (true)
    {
        chrono::time_point current_time = chrono::high_resolution_clock::now();
        Model::dynamic_matrix_t u;
        chrono::duration<double> elapsed = current_time - m_start_time;
        getInterpolatedInput(u, U, elapsed.count() / t);
        commandInput(u);
        std::this_thread::sleep_for(chrono::milliseconds(1));
    }
}
