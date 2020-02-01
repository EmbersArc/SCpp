cmake_minimum_required(VERSION 3.5)
project(scpp_mpc)

# OPTIONS
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(ADDITIONAL_CLEAN_FILES "")

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
file(MAKE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/CppAD/build)
execute_process(COMMAND cmake .. WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/CppAD/build)
list(APPEND ADDITIONAL_CLEAN_FILES ${CMAKE_CURRENT_LIST_DIR}/lib/CppAD/build)

# CppADCodegen
message("\n### Configuring CppADCodegen ###")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/CppADCodegen/build)
execute_process(COMMAND cmake -DCPPAD_INCLUDE_DIRS=${CMAKE_CURRENT_LIST_DIR}/lib/CppAD ..
                WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/CppADCodegen/build)
list(APPEND ADDITIONAL_CLEAN_FILES ${CMAKE_CURRENT_LIST_DIR}/lib/CppADCodegen/build)

# INCLUDES
set(SCPP_SYSTEM_INCLUDE_DIRECTORIES
    lib/external
    lib/CppAD/include
    lib/CppADCodegen/include
    lib/CppADCodegen/build/include
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
    scpp_core/utils/src/parameterServer.cpp
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

add_library(mpc_lib SHARED ${MPC_SOURCES})

target_link_libraries(mpc_lib socp_lib fmt stdc++fs dl Eigen3::Eigen)
target_include_directories(mpc_lib PUBLIC ${SCPP_INCLUDE_DIRECTORIES})
target_include_directories(mpc_lib SYSTEM PUBLIC ${SCPP_SYSTEM_INCLUDE_DIRECTORIES})

set(DEBUG_OPTIONS -Wall -Wextra -Wpedantic)
set(RELEASE_OPTIONS -O2)
target_compile_options(mpc_lib PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")
target_compile_options(mpc_lib PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")


add_executable(mpc_sim scpp/src/MPC_sim.cpp)
target_link_libraries(mpc_sim mpc_lib)

add_executable(lqr_sim scpp/src/LQR_sim.cpp)
target_link_libraries(lqr_sim mpc_lib)

add_executable(sc_tracking scpp/src/SC_tracking.cpp)
target_link_libraries(sc_tracking mpc_lib)

add_executable(sc_oneshot scpp/src/SC_oneshot.cpp)
target_link_libraries(sc_oneshot mpc_lib)
