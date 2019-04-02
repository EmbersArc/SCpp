cmake_minimum_required(VERSION 3.5)
project(SCpp)

# OPTIONS
set(CMAKE_CXX_STANDARD 17)

# ECOS
execute_process(COMMAND make shared WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib/ecos")
include_directories(lib/ecos/include)
include_directories(lib/ecos/external/SuiteSparse_config)
link_directories(lib/ecos)
link_libraries(ecos)

# BOOST
find_package(Boost REQUIRED) 

# FILESYSTEM
link_libraries(stdc++fs)

# FMT
add_subdirectory(lib/fmt EXCLUDE_FROM_ALL)

# CPPAD
include_directories(SYSTEM lib/external)

# KSP
# link_libraries(krpc)
# link_libraries(:libprotobuf.so.17)

# THREAD
link_libraries(pthread)

# INCLUDES
include_directories(models/include)
include_directories(SCpp/include)
include_directories(solver/include)

# SOURCES
set(SOURCES_SC
    models/src/cartpole.cpp
    models/src/rocketLanding3d.cpp
    models/src/rocketLanding2d.cpp
    models/src/crewDragon.cpp
    
    solver/src/ecosWrapper.cpp
    solver/src/optimizationProblem.cpp

    SCpp/src/parameterServer.cpp
    SCpp/src/discretization.cpp
    SCpp/src/successiveConvexificationProblem.cpp
    SCpp/src/timing.cpp

    SCpp/src/freeFinalTimeAlgorithm.cpp
)

add_library(SC STATIC ${SOURCES_SC})
target_compile_options(SC PUBLIC -O2 -Wall -Wextra -Wpedantic -Wno-unused-parameter)

add_executable(SCpp SCpp/src/main.cpp)
# add_executable(SCppKSP SCpp/src/mainKsp.cpp SCpp/src/rocketGuidance.cpp)

target_link_libraries(SC fmt)

target_link_libraries(SCpp SC)
# target_link_libraries(SCppKSP SC)

