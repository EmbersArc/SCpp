cmake_minimum_required(VERSION 3.5)
project(SCpp)

# OPTIONS
set(CMAKE_CXX_STANDARD 17)

# ECOS
execute_process(COMMAND make shared WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib/ecos")
include_directories(socp_interface/solvers/ecos/include)
include_directories(socp_interface/solvers/ecos/external/SuiteSparse_config)
link_directories(socp_interface/solvers/ecos)
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
include_directories(SCpp/include)
include_directories(MPCpp/models/include)
include_directories(SCpp/models/include)
include_directories(socp_interface/include)

# SOURCES
# common utilities
set(UTILS_SOURCES
    utils/src/timing.cpp
    utils/src/parameterServer.cpp
    utils/src/discretization.cpp
    utils/src/simulation
)
# second order cone programming
set(SCOP_SOURCES
    socp_interface/src/ecosWrapper.cpp
    socp_interface/src/optimizationProblem.cpp
)
# MPC algorithm
set(MPC_SOURCES
    MPCpp/src/MPCProblem.cpp
    MPCpp/src/MPCAlgorithm.cpp
)
# SCvx algorithm
set(SCVX_SOURCES
    SCpp/src/SCProblem.cpp
    SCpp/src/SCAlgorithm.cpp
)
# MPC models
set(MPC_MODELS_SOURCES
MPCpp/models/src/rocketHover.cpp
)
# SCvx models
set(SCVX_MODELS_SOURCES
    SCpp/models/src/cartpole.cpp
    SCpp/models/src/rocketLanding3d.cpp
    SCpp/models/src/rocketLanding2d.cpp
    SCpp/models/src/crewDragon.cpp
    SCpp/models/src/rocket3d.cpp
)

add_library(utils_lib STATIC ${UTILS_SOURCES})
add_library(scop_lib STATIC ${SCOP_SOURCES})
add_library(mpc_lib STATIC ${MPC_SOURCES})
add_library(scvx_lib STATIC ${SCVX_SOURCES})
add_library(mpc_models_lib STATIC ${MPC_MODELS_SOURCES})
add_library(scvx_models_lib STATIC ${SCVX_MODELS_SOURCES})

target_link_libraries(utils_lib fmt)
target_link_libraries(scop_lib fmt)
target_link_libraries(mpc_lib scop_lib fmt)
target_link_libraries(scvx_lib scop_lib fmt)
target_link_libraries(mpc_models_lib fmt)
target_link_libraries(scvx_models_lib fmt)

add_executable(MPCpp MPCpp/src/mainMPC.cpp)
add_executable(SCpp SCpp/src/mainSC.cpp)

target_link_libraries(MPCpp scop_lib mpc_lib mpc_models_lib utils_lib)
target_link_libraries(SCpp scop_lib scvx_lib scvx_models_lib utils_lib)
