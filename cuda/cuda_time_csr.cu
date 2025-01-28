#include "csr.h"
#include "utils.h"
#include "vec.h"
#include <stdio.h>
#include <cuda_runtime.h>
#include "cuda_time.h"
#include "spmv_cuda.h"
#include "cuda_mtx.h"

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
    int m = sm.num_rows;
    double *v = d_random(m);
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
    float sum = 0.0;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    for (int i = 0; i < runs_num; i++) {
        cudaEventRecord(start);
        cuda_spmv_csr<<<2, 1024>>>(d_result, d_row_pointer, d_data, d_col_index, d_v, m);
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

