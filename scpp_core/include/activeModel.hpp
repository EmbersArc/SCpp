#pragma once

// --- SC ---

// #include "starship.hpp"
// using Model = scpp::models::Starship;

// --- MPC ---

#include "rocketHover.hpp"
using Model = scpp::models::RocketHover;


using DiscretizationData = Model::DiscretizationData;
using TrajectoryData = Model::TrajectoryData;