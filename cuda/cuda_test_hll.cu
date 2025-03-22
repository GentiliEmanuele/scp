#include "cuda_test.h"
#include "spmv_cuda.h"
#include "cuda_mtx.h"
#include "hll.h"
#include "spmv_openmp.h"
#include "utils.h"
#include "vec.h"
#include <cuda_runtime.h>

//#define hll_v0
//#define hll_v1
#define hll_v2
//#define hll_v3
 

int hll_test(char *path, int hack_size) {
    struct MatrixMarket mm;
    if (read_mtx(path, &mm)) {
        return -1;
    }
    struct hll sm;
    if (hll_init(&sm, hack_size, &mm)) {
        mtx_cleanup(&mm);
        return -1;
    }
    mtx_cleanup(&mm);
    double *v = d_random(mm.num_cols);
    double *d_data;
    int *d_col_index;
    int *d_maxnzr;
    int *d_offsets;
    cudaError_t err = cuda_hll_init(&sm, &d_data, &d_col_index, &d_maxnzr, &d_offsets);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        hll_cleanup(&sm);
        return -1;
    }
    double *d_result;
    err = cudaMalloc(&d_result, sm.num_rows * sizeof(double));
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_maxnzr);
    	hll_cleanup(&sm);
        return 1;
    }
    double *d_v;
    err = cudaMalloc(&d_v, sm.num_rows * sizeof(double));
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_maxnzr);
    	hll_cleanup(&sm);
    	cudaFree(d_result);
        return 1;        
    }
    err = cudaMemcpy(d_v, v, sm.num_rows * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_maxnzr);
    	hll_cleanup(&sm);
    	cudaFree(d_result);
        cudaFree(d_v);
        return 1;
    }
    int threads_num = 1024;
    int blocks_num = (int)ceil(sm.num_rows * 32 / (double)threads_num);
    #ifdef hll_v0
    cuda_spmv_hll_v0<<<blocks_num, threads_num>>>(d_result, sm.hack_size, sm.hacks_num, d_data, d_offsets, d_col_index, d_maxnzr, d_v, sm.num_rows);
    #endif
    #ifdef hll_v1
    cuda_spmv_hll_v1<<<blocks_num, threads_num>>>(d_result, sm.hack_size, sm.hacks_num, d_data, d_offsets, d_col_index, d_maxnzr, d_v, sm.num_rows);
    #endif
    #ifdef hll_v2
    int shared_mem_size = threads_num * sizeof(double);
    cuda_spmv_hll_v2<<<blocks_num, threads_num, shared_mem_size>>>(d_result, sm.hack_size, sm.hacks_num, d_data, d_offsets, d_col_index, d_maxnzr, d_v, sm.num_rows);
    #endif
    #ifdef hll_v3
    cuda_spmv_hll_v3<<<blocks_num, threads_num>>>(d_result, sm.hack_size, sm.hacks_num, d_data, d_offsets, d_col_index, d_maxnzr, d_v, sm.num_rows);
    #endif
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
        cudaFree(d_col_index);
        cudaFree(d_maxnzr);
        cudaFree(d_result);
        cudaFree(d_v);
        hll_cleanup(&sm);
        return -1;
    }
    double *result = d_zeros(sm.num_rows);
    if (result == NULL) {
        printf("cannot allocate result for cuda\n");
    }
    err = cudaMemcpy(result, d_result, sm.num_rows * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
    }
    double *test_result = d_zeros(sm.num_rows);
    if (spmv_hll_par(test_result, &sm, v, sm.num_rows)) {
        printf("spmv_hll_par failed\n");
    } else if (!d_veceq(result, test_result, sm.num_rows, 1e-6)) {
        printf("matrix %s\n", path);
        printf("kernel(%d, %d)\n", blocks_num, threads_num);
        printf("test failed\n");
    } else {
	printf("Test passed for %s \n", path);
    }
    free(test_result);
    free(result);
    cudaFree(d_data);
    cudaFree(d_col_index);
    cudaFree(d_maxnzr);
    cudaFree(d_result);
    cudaFree(d_v);
    free(v);
    hll_cleanup(&sm);
    return 0;
}
