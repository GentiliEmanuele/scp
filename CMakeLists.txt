cmake_minimum_required(VERSION 3.1)
set(CMAKE_CUDA_ARCHITECTURES "native")
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
project(scp LANGUAGES C CXX CUDA)


set(CMAKE_CXX_FLAGS "-O3 -fopenmp -lm")
set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-DSM_${CC};-arch=sm_${CC};-lineinfo;-Xcompiler=${CMAKE_CXX_FLAGS};-Xptxas=-v)

find_package(CUDA REQUIRED)

add_executable(gpu_info cuda/gpuinfo.cu)
add_executable(get_time cuda/cuda_time.cu)
add_executable(do_test cuda/cuda_test.cu)
add_library(time_lib SHARED 
    omp/spmv_openmp.h omp/spmv_openmp.c
    core/csr.h core/csr.c
    core/hll.h core/hll.c
    core/mmio.h core/mmio.c
    core/utils.h core/utils.c
    core/vec.h core/vec.c
    cuda/cuda_time.h cuda/spmv_csr_cuda.cu cuda/cuda_time_hll.cu
    cuda/spmv_cuda.h cuda/spmv_cuda.cu)
target_include_directories(time_lib PUBLIC core cuda omp)
target_link_libraries(get_time PUBLIC time_lib)
target_link_libraries(do_test PUBLIC time_lib)
target_link_libraries(do_test PUBLIC ${CUDA_LIBRARIES})
target_link_libraries(get_time PUBLIC ${CUDA_LIBRARIES})
target_link_libraries(gpu_info PUBLIC ${CUDA_LIBRARIES})
set_target_properties(time_lib PROPERTIES LINKER_LANGUAGE C)

