cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "-O3")
set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-DSM_${CC};-arch=sm_${CC};-lineinfo;-Xcompiler=${CMAKE_CXX_FLAGS};-Xptxas=-v)

find_package(CUDA REQUIRED)

add_executable(test cuda/test.c core/vec.c core/vec.h)
add_library(one STATIC core/vec.h core/vec.c)
target_include_directories(one PUBLIC core)
target_link_libraries(test PUBLIC one)
target_link_libraries(test PUBLIC ${CUDA_LIBRARIES})