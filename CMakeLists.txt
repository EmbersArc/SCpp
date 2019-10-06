cmake_minimum_required(VERSION 3.5)
project(socp_mpc)

# OPTIONS
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Eigen
find_package(Eigen3 REQUIRED NO_MODULE)

# BOOST
find_package(Boost REQUIRED)

# FMT
add_subdirectory(lib/fmt EXCLUDE_FROM_ALL)
set_property(TARGET fmt PROPERTY POSITION_INDEPENDENT_CODE ON)

# SOCP INTERFACE
add_subdirectory(socp_mpc/socp_interface)

# INCLUDES
set(SOCP_SYSTEM_INCLUDE_DIRECTORIES
    lib/external
)
set(SOCP_INCLUDE_DIRECTORIES
    socp_mpc/include
    socp_mpc/models/include
    socp_mpc/utils/include
)

# SOURCES
# MPC Algorithm
set(MPC_SOURCES
    # utilities
    socp_mpc/utils/src/timing.cpp
    socp_mpc/utils/src/parameterServer.cpp

    # main functionality
    socp_mpc/src/simulation.cpp
    socp_mpc/src/discretization.cpp
    socp_mpc/src/MPCProblem.cpp
    socp_mpc/src/MPCAlgorithm.cpp
    socp_mpc/src/SCProblem.cpp
    socp_mpc/src/SCAlgorithm.cpp
    socp_mpc/src/LQR.cpp
    socp_mpc/src/LQRAlgorithm.cpp

    # models
    socp_mpc/models/src/rocketHover.cpp
    socp_mpc/models/src/starship.cpp
)

add_library(mpc_lib SHARED ${MPC_SOURCES})

target_link_libraries(mpc_lib scop_lib fmt stdc++fs)
target_include_directories(mpc_lib PUBLIC ${SOCP_INCLUDE_DIRECTORIES})
target_include_directories(mpc_lib SYSTEM PUBLIC ${SOCP_SYSTEM_INCLUDE_DIRECTORIES})

target_compile_options(mpc_lib PUBLIC -std=c++17 -Wall -Wextra -Wpedantic -O2 -march=native)

add_executable(mpc_sim MPC_sim.cpp)
target_link_libraries(mpc_sim mpc_lib)

add_executable(lqr_sim LQR_sim.cpp)
target_link_libraries(lqr_sim mpc_lib)

add_executable(sc_sim SC_sim.cpp)
target_link_libraries(sc_sim mpc_lib)

add_executable(sc_oneshot SC_oneshot.cpp)
target_link_libraries(sc_oneshot mpc_lib)
