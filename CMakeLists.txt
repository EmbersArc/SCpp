cmake_minimum_required(VERSION 3.5)
project(SCpp)

# OPTIONS
set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter -Werror=return-type -O2 -std=c++17)

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
link_libraries(dl)

# INCLUDES
include_directories(models/include)
include_directories(SCpp/include)
include_directories(solver/include)

# SOURCES
set(SOURCES
    models/src/rocketLanding3d.cpp
    # models/src/rocketLanding2d.cpp
    # models/src/crewDragon.cpp
    
    solver/src/ecosWrapper.cpp
    solver/src/optimizationProblem.cpp

    SCpp/src/parameterServer.cpp
    SCpp/src/discretization.cpp
    SCpp/src/successiveConvexificationSOCP.cpp
    SCpp/src/timing.cpp

    SCpp/src/main.cpp
)

add_executable(SCpp ${SOURCES})
target_link_libraries(SCpp PRIVATE fmt)
