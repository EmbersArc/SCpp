cmake_minimum_required(VERSION 3.5)
project(scpp_mpc)

# OPTIONS
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Eigen
find_package(Eigen3 REQUIRED NO_MODULE)

# BOOST
find_package(Boost REQUIRED)

# FMT
message("\n### Configuring FMT ###")
add_subdirectory(lib/fmt EXCLUDE_FROM_ALL)
set_property(TARGET fmt PROPERTY POSITION_INDEPENDENT_CODE ON)

# SOCP INTERFACE
message("\n### Configuring socp_interface ###")
add_subdirectory(lib/socp_interface)

# CppAD
message("\n### Configuring CppAD ###")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/CppAD)
execute_process(COMMAND cmake ${CMAKE_CURRENT_LIST_DIR}/lib/CppAD 
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/CppAD)

# CppADCodegen
message("\n### Configuring CppADCodegen ###")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/CppADCodegen)
execute_process(COMMAND cmake -DCPPAD_HOME=${CMAKE_CURRENT_LIST_DIR}/lib/CppAD/include ${CMAKE_CURRENT_LIST_DIR}/lib/CppADCodegen
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/CppADCodegen)

# INCLUDES
set(SCPP_SYSTEM_INCLUDE_DIRECTORIES
    lib/external
    lib/CppAD/include
    lib/CppADCodegen/include
    build/lib/CppADCodegen/include
)
set(SCPP_INCLUDE_DIRECTORIES
    scpp_core/include
    scpp_core/utils/include
    scpp_models/include
    scpp/include
)

# SOURCES
# MPC Algorithm
set(MPC_SOURCES
    # utilities
    scpp_core/utils/src/timing.cpp
    scpp/src/commonFunctions.cpp

    # main functionality
    scpp_core/src/simulation.cpp
    scpp_core/src/discretization.cpp
    scpp_core/src/MPCProblem.cpp
    scpp_core/src/MPCAlgorithm.cpp
    scpp_core/src/SCProblem.cpp
    scpp_core/src/SCAlgorithm.cpp
    scpp_core/src/SCvxProblem.cpp
    scpp_core/src/SCvxAlgorithm.cpp
    scpp_core/src/LQR.cpp
    scpp_core/src/LQRAlgorithm.cpp
    scpp_core/src/LQRTracker.cpp

    # models
    scpp_models/src/rocketHover.cpp
    scpp_models/src/rocket2d.cpp
    scpp_models/src/starship.cpp
)

add_library(scpp SHARED ${MPC_SOURCES})

target_link_libraries(scpp socp_interface fmt stdc++fs dl Eigen3::Eigen)
target_include_directories(scpp PUBLIC ${SCPP_INCLUDE_DIRECTORIES})
target_include_directories(scpp SYSTEM PUBLIC ${SCPP_SYSTEM_INCLUDE_DIRECTORIES})

set(DEBUG_OPTIONS -Wall -Wextra -Wpedantic)
set(RELEASE_OPTIONS -O2)
target_compile_options(scpp PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")
target_compile_options(scpp PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")


add_executable(mpc_sim scpp/src/MPC_sim.cpp)
target_link_libraries(mpc_sim scpp)

# add_executable(lqr_sim scpp/src/LQR_sim.cpp)
# target_link_libraries(lqr_sim scpp)

add_executable(sc_tracking scpp/src/SC_tracking.cpp)
target_link_libraries(sc_tracking scpp)

add_executable(sc_oneshot scpp/src/SC_oneshot.cpp)
target_link_libraries(sc_oneshot scpp)
