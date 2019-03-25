#ifndef ROCKET_GUIDANCE_HPP
#define ROCKET_GUIDANCE_HPP

#include <chrono>
#include <thread>

#include <krpc.hpp>
#include <krpc/services/krpc.hpp>
#include <krpc/services/space_center.hpp>
#include <krpc/services/ui.hpp>
#include <krpc/services/drawing.hpp>

#include "freeFinalTimeAlgorithm.hpp"
#include "timing.hpp"

using namespace krpc::services;

typedef std::tuple<double, double, double> tuple3d_t;
typedef std::tuple<double, double, double, double> tuple4d_t;

class RocketGuidance
{
  public:
    RocketGuidance();
    void run();

  private:
    void initializeReferenceFrames();
    void drawAxes();
    void drawTrajectory(const Eigen::MatrixXd &X);
    void prepareVessel();
    void updateModelParameters();
    void recomputeTrajectory(Eigen::MatrixXd &X, Eigen::MatrixXd &U, double &t, bool &running);
    void commandInput(const Eigen::MatrixXd &u);
    static void getInterpolatedInput(Eigen::MatrixXd &u, const Eigen::MatrixXd &U, const double t);

    // services
    krpc::Client m_connection;
    std::shared_ptr<KRPC> m_krpc;
    std::shared_ptr<SpaceCenter> m_spaceCenter;
    std::shared_ptr<Drawing> m_drawing;

    // vessel and celestial body
    SpaceCenter::Vessel m_vessel;
    SpaceCenter::CelestialBody m_body;

    // reference frames
    SpaceCenter::ReferenceFrame m_ref_frame_pad;
    SpaceCenter::ReferenceFrame m_ref_frame_vessel;
    SpaceCenter::ReferenceFrame m_ref_frame_nonrot;

    // streams
    krpc::Stream<float> m_m_stream;
    krpc::Stream<tuple3d_t> m_r_I_stream;
    krpc::Stream<tuple3d_t> m_v_I_stream;
    krpc::Stream<tuple4d_t> m_q_I_stream;
    krpc::Stream<tuple3d_t> m_w_B_stream;
    krpc::Stream<tuple3d_t> m_J_B_stream;
    krpc::Stream<tuple3d_t> m_r_T_B_stream;
    krpc::Stream<double> m_t_stream;
    krpc::Stream<float> m_pitch_stream;
    krpc::Stream<float> m_yaw_stream;
    krpc::Stream<float> m_T_stream;

    // engine reference
    SpaceCenter::Engine m_engine;

    // model and algorithm
    std::shared_ptr<Model> m_model;
    std::shared_ptr<FreeFinalTimeAlgorithm> m_solver;
    std::atomic<bool> m_warm_start;
    std::thread m_solverThread;
    std::thread m_drawingThread;
    std::chrono::system_clock::time_point m_start_time;
};

#endif