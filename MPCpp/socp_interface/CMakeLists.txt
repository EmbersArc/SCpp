cmake_minimum_required(VERSION 3.5)
project(socp_interface)

# OPTIONS
set(CMAKE_CXX_STANDARD 17)

# ECOS
execute_process(COMMAND make shared WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/solvers/ecos")
# include_directories(solvers/ecos/include)
# include_directories(solvers/ecos/external/SuiteSparse_config)
# link_directories(solvers/ecos)

set(SCOP_SOURCES
    src/ecosWrapper.cpp
    src/optimizationProblem.cpp
)

add_library(socp_interface_lib STATIC ${SCOP_SOURCES})
target_include_directories(socp_interface_lib PUBLIC ${CMAKE_CURRENT_LIST_DIR}/solvers/ecos/include ${CMAKE_CURRENT_LIST_DIR}/solvers/ecos/external/SuiteSparse_config ${CMAKE_CURRENT_LIST_DIR}/include)
