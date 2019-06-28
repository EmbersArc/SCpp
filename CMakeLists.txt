cmake_minimum_required(VERSION 3.5)
project(socp_mpc)

# OPTIONS
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -O2")

# Eigen
find_package(Eigen3 REQUIRED NO_MODULE)

# BOOST
find_package(Boost REQUIRED)

# FMT
add_subdirectory(lib/fmt EXCLUDE_FROM_ALL)

# CPPAD
include_directories(SYSTEM lib/external)

# SOCP INTERFACE
add_subdirectory(socp_mpc/socp_interface)

# INCLUDES
include_directories(socp_mpc/include)
include_directories(socp_mpc/models/include)
include_directories(socp_mpc/utils/include)

# SOURCES
# MPC Algorithm
set(MPC_SOURCES
    socp_mpc/utils/src/timing.cpp
    socp_mpc/utils/src/parameterServer.cpp

    socp_mpc/src/simulation.cpp
    socp_mpc/src/discretization.cpp
    socp_mpc/src/MPCProblem.cpp
    socp_mpc/src/MPCAlgorithm.cpp
    socp_mpc/src/SCProblem.cpp
    socp_mpc/src/SCAlgorithm.cpp

    socp_mpc/models/src/rocketHover.cpp
    socp_mpc/models/src/rocketLanding3d.cpp
)

add_library(mpc_lib STATIC ${MPC_SOURCES})

target_link_libraries(mpc_lib scop_lib fmt)

add_executable(mpc mainMPC.cpp)
target_link_libraries(mpc mpc_lib)
add_executable(sc mainSC.cpp)
target_link_libraries(sc mpc_lib)
