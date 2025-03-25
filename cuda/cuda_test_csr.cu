#include "cuda_test.h"
#include "spmv_cuda.h"
#include "cuda_mtx.h"
#include "csr.h"
#include "spmv_openmp.h"
#include "utils.h"
#include "vec.h"
#include <cuda_runtime.h>
#include <math.h>

#define WARP_SIZE 32

//#define csr_v1
//#define csr_v2
//#define csr_v3
#define csr_v4

int csr_test(const char *path) {
    struct MatrixMarket mm;
    if (read_mtx(path, &mm)) {
        return -1;
    }
    struct csr sm;
    if (csr_init(&sm, &mm)) {
        mtx_cleanup(&mm);
        return -1;
    }
    mtx_cleanup(&mm);
    double *d_data;
    int *d_col_index;
    int *d_row_pointer;
    cudaError_t err = cuda_csr_init(&sm, &d_data, &d_col_index, &d_row_pointer);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        csr_cleanup(&sm);
        return -1;
    }
    double *d_result;
    err = cudaMalloc(&d_result, sm.num_rows * sizeof(double));
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
        cudaFree(d_row_pointer);
    	csr_cleanup(&sm);
        return 1;
    }
    double *d_v;
    err = cudaMalloc(&d_v, sm.num_rows * sizeof(double));
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
        cudaFree(d_row_pointer);
    	cudaFree(d_result);
    	csr_cleanup(&sm);
        return 1;        
    }
    double *v = d_random(sm.num_cols);
    err = cudaMemcpy(d_v, v, sm.num_rows * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_row_pointer);
    	cudaFree(d_result);
        cudaFree(d_v);
        csr_cleanup(&sm);
        return 1;
    }
    #ifdef csr_v1
    int threads_num = 1024;
    int blocks_num = (int)ceil(sm.num_rows / (double)threads_num);
    cuda_spmv_csr<<<blocks_num, threads_num>>>(d_result, d_row_pointer, d_data, d_col_index, d_v, sm.num_rows);
    #endif
    #ifdef csr_v2
    int threads_num = 1024;
    int blocks_num = (int)ceil(sm.num_rows / (double)threads_num);
    int shared_mem_size = threads_num * sizeof(double);
    cuda_spmv_csr_v2<<<blocks_num, threads_num, shared_mem_size>>>(d_result, d_row_pointer, d_data, d_col_index, d_v, sm.num_rows);
    #endif
    #ifdef csr_v3
    int threads_num = 1024;
    int blocks_num = (int)ceil(sm.num_rows *(double) WARP_SIZE / (double)threads_num);
    cuda_spmv_csr_v3<<<blocks_num, threads_num>>>(d_result, d_row_pointer, d_data, d_col_index, d_v, sm.num_rows);
    #endif
    #ifdef csr_v4
    int threads_num = 1024;
    int blocks_num = (int)ceil(sm.num_rows *(double) WARP_SIZE / (double)threads_num);
    int shared_mem_size = threads_num * sizeof(double);
    cuda_spmv_csr_v4<<<blocks_num, threads_num, shared_mem_size>>>(d_result, d_row_pointer, d_data, d_col_index, d_v, sm.num_rows);
    #endif
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
        cudaFree(d_col_index);
        cudaFree(d_result);
        cudaFree(d_row_pointer);
        cudaFree(d_v);
        csr_cleanup(&sm);
        return 1;
    }
    double *result = d_zeros(sm.num_rows);
    err = cudaMemcpy(result, d_result, sm.num_rows * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
    }
    double *test_result = d_zeros(sm.num_rows);
    if (spmv_csr_par(test_result, &sm, v, sm.num_rows)) {
        printf("spmv_csr_par failed\n");
    } else if (!d_veceq(result, test_result, sm.num_rows, 1e-6)) {
        printf("Test failed for %s (#threads=%d, #rows=%d)\n", path, threads_num * blocks_num, sm.num_rows);
    } else {
	printf("Test passed for %s \n", path);
    }
    cudaFree(d_data);
    cudaFree(d_col_index);
    cudaFree(d_result);
    cudaFree(d_row_pointer);
    cudaFree(d_v);
    csr_cleanup(&sm);
    free(test_result);
    free(result);
    free(v);
    return 0;
}
