cmake_minimum_required(VERSION 3.5)
project(SCpp)

# OPTIONS
set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter -Werror=return-type -O1 -std=c++17)

# ECOS
execute_process(COMMAND make shared WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib/ecos")
include_directories(lib/ecos/include)
include_directories(lib/ecos/external/SuiteSparse_config)
link_directories(lib/ecos)
link_libraries(ecos)

# BOOST
find_package(Boost REQUIRED) 

# OMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# FILESYSTEM
link_libraries(stdc++fs)

# FMT
add_subdirectory(lib/fmt EXCLUDE_FROM_ALL)

# CPPAD
link_libraries(dl)

# INCLUDES
include_directories(include)

# SOURCES
set(SOURCES
    include/models/modelRocketLanding3dDefinitions.hpp
    include/models/modelRocketLanding3d.hpp

    include/ecosWrapper.hpp
        src/ecosWrapper.cpp
    include/eigenIntegration.hpp

    include/systemModel.hpp
    include/optimizationProblem.hpp
        src/optimizationProblem.cpp
    include/discretization.hpp
        src/discretization.cpp
    include/successiveConvexificationSOCP.hpp
        src/successiveConvexificationSOCP.cpp
    include/timing.hpp
        src/timing.cpp

        src/main.cpp
)

add_executable(SCpp ${SOURCES})
target_link_libraries(SCpp PRIVATE fmt-header-only)

