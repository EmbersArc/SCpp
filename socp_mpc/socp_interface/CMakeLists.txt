cmake_minimum_required(VERSION 3.5)
project(scop_interface)

# OPTIONS
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -O2")

# ECOS
add_custom_target(ecos 
                  COMMAND make shared
                  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/solvers/ecos)
link_directories(solvers/ecos)
link_libraries(${CMAKE_CURRENT_LIST_DIR}/solvers/ecos/libecos.so)

set(TARGET_INCLUDE
solvers/ecos/include
solvers/ecos/external/SuiteSparse_config
include
)

# SOCP Interface
set(SCOP_SOURCES
    src/ecosWrapper.cpp
    src/optimizationProblem.cpp
)

add_library(scop_lib SHARED ${SCOP_SOURCES})
target_include_directories(scop_lib PUBLIC ${TARGET_INCLUDE})