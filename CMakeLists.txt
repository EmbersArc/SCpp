cmake_minimum_required(VERSION 3.9)
project(OptimalLanding)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-Wall -fopenmp -O3 " )

add_executable(OptimalLanding
    include/params.h
    src/main.cpp)

include_directories(./include)

add_compile_options()