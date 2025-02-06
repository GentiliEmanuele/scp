#include "csr.h"
#include "utils.h"
#include "vec.h"
#include "cuda_time.h"
#include "spmv_cuda.h"
#include "cuda_mtx.h"
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>

int csr_time(const char *path, int runs_num, struct time_info *ti) {
    struct MatrixMarket mm;
    if (read_mtx(path, &mm)) {
        return -1;
    }
    struct csr sm;
    if (csr_init(&sm, &mm)) {
        mtx_cleanup(&mm);
        return -1;
    }
    int nz = mm.nz;
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
#ifdef CUDA_MANAGED
    err = cudaMallocManaged(&d_result, sm.num_rows * sizeof(double));
#else
    err = cudaMalloc(&d_result, sm.num_rows * sizeof(double));
#endif
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
        cudaFree(d_row_pointer);
    	csr_cleanup(&sm);
        return 1;
    }
    double *d_v;
#ifdef CUDA_MANAGED
    err = cudaMallocManaged(&d_v, sm.num_cols * sizeof(double));
#else
    err = cudaMalloc(&d_v, sm.num_cols * sizeof(double));
#endif
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
    err = cudaMemcpy(d_v, v, sm.num_cols * sizeof(double), cudaMemcpyHostToDevice);
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
    float sum = 0.0;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    int threads_num = 1024;
    int blocks_num = (int)ceil(sm.num_rows / (double)threads_num);
    for (int i = 0; i < runs_num; i++) {
        cudaEventRecord(start);
        cuda_spmv_csr<<<blocks_num, threads_num>>>(d_result, d_row_pointer, d_data, d_col_index, d_v, sm.num_rows);
        cudaEventRecord(stop);
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
        cudaEventSynchronize(stop);
        float m = 0.0;
        cudaEventElapsedTime(&m, start, stop);
        sum += m;
    }
    double *result = d_zeros(sm.num_rows);
    err = cudaMemcpy(result, d_result, sm.num_rows * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
    }
    ti->millis = sum / runs_num;
    ti->flops = (2 * nz) / ti->millis;
    cudaFree(d_data);
    cudaFree(d_col_index);
    cudaFree(d_result);
    cudaFree(d_row_pointer);
    cudaFree(d_v);
    csr_cleanup(&sm);
    free(v);
    return 0;
}

