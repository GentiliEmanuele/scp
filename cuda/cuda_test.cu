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

void print_vec(double *v, int n) {
    for (int i = 0; i < n; i++)
    {
        printf("%d %f\n", i, v[i]);
    }
}

int parse_test_type(char *s) {
    if (!strcmp(s, "csr")) {
        return CSR;
    } else if (!strcmp(s, "hll")) {
        return HLL;
    } else {
        return NOTEST;
    }
}

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
    } else if (d_veceq(result, test_result, sm.num_rows, 1e-6)) {
        printf("test passed\n");
    } else {
        printf("test failed\n");
    }
    printf("cuda result\n");
    print_vec(result, 10);
    printf("oracle result\n");
    print_vec(test_result, 10);
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

    int test_type = parse_test_type(argv[2]);
    if (test_type == NOTEST) {
        printf("expected one of csr or hll but got %s\n", argv[2]);
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
