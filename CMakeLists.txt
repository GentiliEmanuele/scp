cmake_minimum_required(VERSION 3.1)
set(CMAKE_CUDA_ARCHITECTURES "native")
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
project(scp LANGUAGES C CXX CUDA)


set(CMAKE_CXX_FLAGS "-O3 -fopenmp -lm")
set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-DSM_${CC};-arch=sm_${CC};-lineinfo;-Xcompiler=${CMAKE_CXX_FLAGS};-Xptxas=-v)

find_package(CUDA REQUIRED)

add_executable(gpu_info cuda/gpuinfo.cu)
add_executable(csr_time cuda/spmv_csr_cuda.cu)
add_executable(hll_time cuda/cuda_time_hll.cu)
add_library(csr_time_lib SHARED core/csr.h core/csr.c core/utils.h core/utils.c core/vec.h core/vec.c cuda/cuda_time.h cuda/spmv_cuda.h)
add_library(hll_time_lib SHARED cuda/cuda_time.h core/hll.h cuda/spmv_cuda.h core/utils.h core/utils.c core/vec.h core/vec.c)
target_include_directories(csr_time_lib core cuda)
target_include_directories(hll_time_lib core cuda)
target_link_libraries(csr_time PUBLIC csr_time_lib)
target_link_libraries(hll_time PUBLIC hll_time_lib)
target_link_libraries(csr_time PUBLIC ${CUDA_LIBRARIES})
target_link_libraries(hll_time PUBLIC ${CUDA_LIBRARIES})
target_link_libraries(gpu_info PUBLIC ${CUDA_LIBRARIES})
set_target_properties(csr_time_lib PROPERTIES LINKER_LANGUAGE C)
set_target_properties(hll_time_lib PROPERTIES LINKER_LANGUAGE C)

