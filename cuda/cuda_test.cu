#include "cuda_time.h"
#include "cuda_mtx.h"
#include "spmv_cuda.h"
#include "spmv_openmp.h"
#include "csr.h"
#include "hll.h"
#include "utils.h"
#include "vec.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define pr_err(err) printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err))

#define NOTEST 0
#define CSR    1
#define HLL    2

int parse_test_type(char *s) {
    if (!strcmp(s, "csr")) {
        return CSR;
    } else if (!strcmp(s, "hll")) {
        return HLL;
    } else {
        return NOTEST;
    }
}

int csr_test(char *path) {
    struct MatrixMarket mm;
    if (read_mtx(path, &mm)) {
        return 1;
    }

    struct csr sm;
    if (csr_init(&sm, &mm)) {
        mtx_cleanup(&mm);
        return 1;
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
    double *par_result = d_zeros(m);
    double *gpu_result = d_zeros(m);
    cudaMemcpy(gpu_result, d_result, sizeof(double) * m, cudaMemcpyDeviceToHost);
    if(spmv_csr_par(par_result, &sm, v, m)) {
        printf("An error occurred executing spmv_csr_par");
    } else if (d_veceq(par_result, gpu_result, sm.num_rows, 1e-6)) {
        printf("test passed\n");
    } else {
        printf("test failed\n");
    }
    
    cudaFree(d_data);
    cudaFree(d_col_index);
    cudaFree(d_result);
    cudaFree(d_row_pointer);
    cudaFree(d_v);
    csr_cleanup(&sm);
    free(v);
    return 0;
}

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
    cuda_spmv_hll<<<1024, 1024>>>(d_result, sm.hack_size, sm.hacks_num, d_data, d_offsets, d_col_index, d_maxnzr, d_v, sm.num_rows);
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
    double *result = d_zeros(sm.num_rows);
    cudaMemcpy(result, d_result, sm.num_rows * sizeof(double), cudaMemcpyDeviceToHost);
    
    double *test_result = d_zeros(sm.num_rows);
    if (spmv_hll_par(test_result, &sm, v, sm.num_rows)) {
        printf("spmv_hll_par failed\n");
    } else if (d_veceq(result, test_result, sm.num_rows, 1e-6)) {
        printf("test passed\n");
    } else {
        printf("test failed\n");
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

int main(int argc, char *argv[]) {
    --argc;
    if (argc != 2 && argc != 3) {
        printf("see usage: program matrix [hll|csr] {hack_size}\n");
        return -1;
    }

    int test_type = parse_test_type(argv[4]);
    if (test_type == NOTEST) {
        printf("expected one of csr or hll but got %s\n", argv[4]);
        return -1;
    }

    int hack_size = 0;
    if (test_type == HLL) {
        if (argc == 3) {
            hack_size = atoi(argv[3]);
            if (hack_size == 0) {
                printf("An error occurred while converting hack_size\n");
                return -1;
            }
        } else {
            hack_size = 32;
            printf("no hack_size specified, default value (32) set\n");
        }
    }
    
    printf("matrix %s\n", argv[1]);
    if (test_type == CSR) {
        csr_test(argv[1]);
    } else if (test_type == HLL) {
        hll_test(argv[1], hack_size);
    }
}
