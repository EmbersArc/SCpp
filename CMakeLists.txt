cmake_minimum_required(VERSION 3.5)
project(OptimalLanding)

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-Wall -Wextra -Wpedantic -O1) 

find_package(Boost REQUIRED)   
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

include_directories(./include)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(OptimalLanding
    include/params.h
    src/main.cpp)
