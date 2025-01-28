#include "cuda_test.h"
#include "cuda_mtx.h"
#include "csr.h"
#include "spmv_openmp.h"
#include "utils.h"
#include "vec.h"
#include <cuda_runtime.h>

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
    double *v = d_random(sm.num_rows);
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
    for (int i = 0; i < 1; i++) {
        cuda_spmv_csr<<<2, 1024>>>(d_result, d_row_pointer, d_data, d_col_index, d_v, m);
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
    }
    double *result = d_zeros(sm.num_rows);
    err = cudaMemcpy(result, d_result, sm.num_rows * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
    }
    print_vec(result, 10);
    cudaFree(d_data);
    cudaFree(d_col_index);
    cudaFree(d_result);
    cudaFree(d_row_pointer);
    cudaFree(d_v);
    csr_cleanup(&sm);
    free(v);
    return 0;
}