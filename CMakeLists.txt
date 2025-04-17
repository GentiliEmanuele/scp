cmake_minimum_required(VERSION 3.1)
set(CMAKE_CUDA_ARCHITECTURES 75)
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
project(scp LANGUAGES C CXX CUDA)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Iinclude -O3 -march=native -fopenmp")
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -lm")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --compiler-options=-fopenmp")


add_compile_options(-ffast-math -funroll-loops -Dcsr_v2 -Dhll_v2)

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
    cuda/cuda_time.h cuda/cuda_time_csr.cu cuda/cuda_time_hll.cu
    cuda/spmv_cuda.h cuda/spmv_cuda.cu
    cuda/cuda_test.h cuda/cuda_test_csr.cu cuda/cuda_test_hll.cu
    cuda/cuda_mtx.h cuda/cuda_mtx.cu)
target_include_directories(time_lib PUBLIC core cuda omp)
target_link_libraries(get_time PUBLIC time_lib)
target_link_libraries(do_test PUBLIC time_lib)
target_link_libraries(do_test PUBLIC ${CUDA_LIBRARIES})
target_link_libraries(get_time PUBLIC ${CUDA_LIBRARIES})
target_link_libraries(gpu_info PUBLIC ${CUDA_LIBRARIES})
set_target_properties(time_lib PROPERTIES LINKER_LANGUAGE C)

