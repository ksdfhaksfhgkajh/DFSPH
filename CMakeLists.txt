cmake_minimum_required(VERSION 3.17)
project(DFSPH LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_ARCHITECTURES 75)

file(GLOB CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB CU_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(solver ${CPP_SOURCES} ${CU_SOURCES}
        src/DFSPHSolver.cpp
)

set_target_properties(solver PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin
)

