#include "cuda_time.h"
#include "cuda_mtx.h"
#include "hll.h"
#include "spmv_cuda.h"
#include "utils.h"
#include "vec.h"
#include <stdio.h>
#include <cuda_runtime.h>

#define pr_err(err) printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err))

int hll_time(const char *path, int runs_num, int hack_size, struct time_info *ti) {
    struct MatrixMarket mm;
    if (read_mtx(path, &mm)) {
        return -1;
    }
    struct hll sm;
    if (hll_init(&sm, hack_size, &mm)) {
        mtx_cleanup(&mm);
        return -1;
    }
    int nz = mm.nz;
    mtx_cleanup(&mm);
    double *v = d_random(mm.num_rows);
    double *d_data;
    int *d_col_index;
    int *d_maxnzr;
    int *d_offsets;
    cudaError_t err = cuda_hll_init(&sm, &d_data, &d_col_index, &d_maxnzr, &d_offsets);
    if (err != cudaSuccess) {
        pr_err(err);
        hll_cleanup(&sm);
        return -1;
    }
    double *d_result;
    err = cudaMalloc(&d_result, sm.num_rows * sizeof(double));
    if (err != cudaSuccess) {
        pr_err(err);
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_maxnzr);
    	hll_cleanup(&sm);
        return 1;
    }
    double *d_v;
    err = cudaMalloc(&d_v, sm.num_rows * sizeof(double));
    if (err != cudaSuccess) {
        pr_err(err);
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_maxnzr);
    	hll_cleanup(&sm);
    	cudaFree(d_result);
        return 1;        
    }
    err = cudaMemcpy(d_v, v, sm.num_rows * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        pr_err(err);
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_maxnzr);
    	hll_cleanup(&sm);
    	cudaFree(d_result);
        cudaFree(d_v);
        return 1;
    }
    float sum = 0.0;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    for (int i = 0; i < runs_num; ++i) {
        cudaEventRecord(start);
        cuda_spmv_hll<<<2, 1024>>>(d_result, sm.hack_size, sm.hacks_num, d_data, d_offsets, d_col_index, d_maxnzr, d_v, sm.num_rows);
        cudaEventRecord(stop);
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            pr_err(err);
            cudaFree(d_data);
            cudaFree(d_col_index);
            cudaFree(d_maxnzr);
            cudaFree(d_result);
            cudaFree(d_v);
            hll_cleanup(&sm);
            return -1;
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
    cudaFree(d_maxnzr);
    cudaFree(d_result);
    cudaFree(d_v);
    free(v);
    hll_cleanup(&sm);
    return 0;
}