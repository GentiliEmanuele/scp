cmake_minimum_required(VERSION 3.1)
set(CMAKE_CUDA_ARCHITECTURES "native")
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
project(scp LANGUAGES C CXX CUDA)


set(CMAKE_CXX_FLAGS "-O3 -fopenmp -lm")
set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-DSM_${CC};-arch=sm_${CC};-lineinfo;-Xcompiler=${CMAKE_CXX_FLAGS};-Xptxas=-v)

find_package(CUDA REQUIRED)

add_executable(gpu_info cuda/gpuinfo.cu)
add_executable(test cuda/spmv_cuda.cu)
add_library(one SHARED core/vec.h core/vec.c core/csr.h core/csr.c core/utils.h core/utils.c core/mmio.h core/mmio.c core/hll.h core/hll.c)
add_library(two SHARED omp/spmv_openmp.h omp/spmv_openmp.c)
target_include_directories(one PUBLIC core)
target_include_directories(two PUBLIC omp)
target_link_libraries(two PUBLIC one)
target_link_libraries(test PUBLIC one)
target_link_libraries(test PUBLIC two)
target_link_libraries(test PUBLIC ${CUDA_LIBRARIES})
target_link_libraries(gpu_info PUBLIC ${CUDA_LIBRARIES})
set_target_properties(one PROPERTIES LINKER_LANGUAGE C)
