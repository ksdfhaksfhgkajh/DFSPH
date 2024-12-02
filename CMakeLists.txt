cmake_minimum_required(VERSION 3.17)
project(DFSPH)

enable_language(CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_ARCHITECTURES 75)

file(GLOB SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu"
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(solver ${SRC})

set_target_properties(solver PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin
)
