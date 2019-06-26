cmake_minimum_required(VERSION 3.5)
project(MPCpp)

# OPTIONS
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2")

# ECOS
add_custom_target(ecos COMMAND make shared
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib/ecos)
include_directories(lib/ecos/include)
include_directories(lib/ecos/external/SuiteSparse_config)
link_directories(lib/ecos)
link_libraries(ecos)

# Eigen
find_package(Eigen3 REQUIRED NO_MODULE)

# BOOST
find_package(Boost REQUIRED)

# FMT
add_subdirectory(lib/fmt EXCLUDE_FROM_ALL)

# CPPAD
include_directories(SYSTEM lib/external)

# INCLUDES
include_directories(utils/include)
include_directories(MPCpp/include)
include_directories(MPCpp/models/include)
include_directories(MPCpp/socp_interface/include)
include_directories(MPCpp/utils/include)

# SOURCES
# SOCP Interface
set(SCOP_SOURCES
    MPCpp/socp_interface/src/ecosWrapper.cpp
    MPCpp/socp_interface/src/optimizationProblem.cpp
)

# MPC Algorithm
set(MPC_SOURCES
    MPCpp/utils/src/timing.cpp
    MPCpp/utils/src/parameterServer.cpp

    MPCpp/src/simulation.cpp
    MPCpp/src/discretization.cpp
    MPCpp/src/MPCProblem.cpp
    MPCpp/src/MPCAlgorithm.cpp

    MPCpp/models/src/rocketHover.cpp
)

add_library(scop_lib STATIC ${SCOP_SOURCES})
add_library(mpc_lib STATIC ${MPC_SOURCES})

target_link_libraries(mpc_lib scop_lib fmt)

add_executable(mpc mainMPC.cpp)
target_link_libraries(mpc mpc_lib)
