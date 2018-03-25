cmake_minimum_required(VERSION 3.5)
project(OptimalLanding)

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-Wall -Wextra -Wpedantic -O1)

find_package(Boost REQUIRED)   
find_package(OpenMP REQUIRED)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

include_directories(./include)
include_directories(eigen)
include_directories(ecos/include)
include_directories(ecos/external/SuiteSparse_config)

link_directories(ecos)
link_libraries(ecos)

add_executable(OptimalLanding
    src/main.cpp
    include/model_landing_6dof.h
    src/model_landing_6dof.cpp
    include/model_simple_4th_order.hpp
    src/model_simple_4th_order.cpp
)
