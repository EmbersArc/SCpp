cmake_minimum_required(VERSION 3.5)
project(OptimalLanding)

# Build ECOS shared library
execute_process(COMMAND make shared WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib/ecos")

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter -Werror=return-type -O1 -std=c++17)

find_package(Boost REQUIRED)   

include_directories(include)
include_directories(lib/ecos/include)
include_directories(lib/ecos/external/SuiteSparse_config)
include_directories(SYSTEM lib/eigen)

link_directories(lib/ecos)
link_libraries(ecos)

set(SOURCES
    include/eigenIntegration.hpp
    include/models/model_landing_6dof.hpp
        src/models/model_landing_6dof.cpp
    # include/models/model_landing_3dof.hpp
    #     src/models/model_landing_3dof.cpp
    # include/models/model_simple_4th_order.hpp
    #     src/models/model_simple_4th_order.cpp

    include/ecosWrapper.hpp
        src/ecosWrapper.cpp
        
    # include/check_jacobian.hpp
    #     src/check_jacobian.cpp
    # include/check_trajectory.hpp
    #     src/check_trajectory.cpp

    include/optimizationProblem.hpp
        src/optimizationProblem.cpp
    include/discretization.hpp
        src/discretization.cpp
    include/successiveConvexificationSOCP.hpp
        src/successiveConvexificationSOCP.cpp
    include/timing.hpp
        src/timing.cpp
)


### Use MOSEK if available
### You need to set the environment variable: 'export MOSEK_ROOT=.../mosek/8/tools/platform/linux64x86'
#if(DEFINED ENV{MOSEK_ROOT})
#    message("MOSEK_ROOT: " $ENV{MOSEK_ROOT})
#    include_directories(SYSTEM $ENV{MOSEK_ROOT}/h)
#    link_directories($ENV{MOSEK_ROOT}/bin)
#    link_libraries(mosek64 fusion64)
#    set(SOURCES ${SOURCES}
#        include/MosekWrapper.hpp
#            src/MosekWrapper.cpp
#    )
#endif()


add_library(successiveConvexification STATIC ${SOURCES})
link_libraries(successiveConvexification)

add_executable(OptimalLanding    src/main.cpp)
# add_executable(test_OptimizationProblem_indices    test/test_OptimizationProblem_indices.cpp)
# add_executable(test_ecos_wrapper    test/test_ecos_wrapper.cpp)
# add_executable(test_ecos_wrapper_sparse_conversion    test/test_ecos_wrapper_sparse_conversion.cpp)
# add_executable(test_check_jacobian    test/test_check_jacobian.cpp)
