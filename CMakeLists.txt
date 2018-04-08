cmake_minimum_required(VERSION 3.5)
project(OptimalLanding)


# Build ECOS shared library
execute_process(COMMAND make shared WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/ecos")


set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter -Werror=return-type -O1)

find_package(Boost REQUIRED)   
find_package(OpenMP REQUIRED)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

include_directories(include)
include_directories(ecos/include)
include_directories(ecos/external/SuiteSparse_config)
include_directories(SYSTEM eigen)

link_directories(ecos)
link_libraries(ecos)

set(SOURCES
    include/model_landing_6dof.h
        src/model_landing_6dof.cpp
    include/model_landing_3dof.h
        src/model_landing_3dof.cpp
    include/model_simple_4th_order.hpp
        src/model_simple_4th_order.cpp
    include/EcosWrapper.hpp
        src/EcosWrapper.cpp
    include/check_jacobian.h
        src/check_jacobian.cpp
    include/check_trajectory.h
        src/check_trajectory.cpp
    include/OptimizationProblem.hpp
        src/OptimizationProblem.cpp
    include/Discretization.hpp
        src/Discretization.cpp
    include/SuccessiveConvexificationSOCP.hpp
        src/SuccessiveConvexificationSOCP.cpp
    include/timing.hpp
        src/timing.cpp
)

## Use MOSEK if available
## You need to set the environment variable: 'export MOSEK_ROOT=.../mosek/8/tools/platform/linux64x86'
if(DEFINED ENV{MOSEK_ROOT})
    message("MOSEK_ROOT: " $ENV{MOSEK_ROOT})
    include_directories(SYSTEM $ENV{MOSEK_ROOT}/h)
    link_directories($ENV{MOSEK_ROOT}/bin)
    link_libraries(mosek64 fusion64)
    set(SOURCES ${SOURCES}
        include/MosekWrapper.hpp
            src/MosekWrapper.cpp
    )
endif()


add_library(successiveConvexification STATIC ${SOURCES})
link_libraries(successiveConvexification)

add_executable(OptimalLanding    src/main.cpp)
add_executable(test_OptimizationProblem_indices    test/test_OptimizationProblem_indices.cpp)
add_executable(test_ecos_wrapper    test/test_ecos_wrapper.cpp)
add_executable(test_ecos_wrapper_sparse_conversion    test/test_ecos_wrapper_sparse_conversion.cpp)
add_executable(test_check_jacobian    test/test_check_jacobian.cpp)
