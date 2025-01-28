#include "cuda_time.h"
#include "cuda_mtx.h"
#include "spmv_cuda.h"
#include "csr.h"
#include "hll.h"
#include "utils.h"
#include "vec.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

void print_vec(double *v, int n) {
    for (int i = 0; i < n; i++)
    {
        printf("%d %f\n", i, v[i]);
    }
}

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
    double *result = d_zeros(sm.num_rows);
    err = cudaMemcpy(result, d_result, sm.num_rows * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
    }
    print_vec(result, 10);
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

int main(int argc, char *argv[]) {
    --argc;
    if (argc != 4 && argc != 5) {
        printf("see usage: program output matrices_list runs_num [hll|csr] {hack_size}\n");
        return -1;
    }

    int runs_num = atoi(argv[3]);
    if (runs_num == 0) {
        printf("An error occurred while converting hack_size\n");
        return -1;
    }
    int test_type = parse_test_type(argv[4]);
    if (test_type == NOTEST) {
        printf("expected one of csr or hll but got %s\n", argv[4]);
        return -1;
    }

    int hack_size = 0;
    if (test_type == HLL) {
        if (argc == 5) {
            hack_size = atoi(argv[5]);
            if (hack_size == 0) {
                printf("An error occurred while converting hack_size\n");
                return -1;
            }
        } else {
            hack_size = 32;
            printf("no hack_size specified, default value (32) set\n");
        }
    }
    
    char filename[1024];
    if (hack_size) {
        sprintf(filename, "%s_%d_%d.csv", argv[1], runs_num, hack_size);
    } else {
        sprintf(filename, "%s_%d.csv", argv[1], runs_num);
    }
    FILE* off = fopen(filename, "w");
    if (off == NULL) {
        printf("cannot open file %s\n", filename);
        return -1;
    }
    FILE* iff = fopen(argv[2], "r");
    if (iff == NULL) {
        printf("cannot open file %s\n", argv[2]);
        fclose(off);
        return -1;
    }
    fprintf(off, "matrix,time,flops\n");
    INIT_TIME_INFO(ti);
    char line[1024];
    while (fgets(line, 1024, iff) != NULL) {
        int n = strlen(line);
        if (line[--n] == '\n') {
            line[n] = 0;
        }
        printf("matrix %s\n", line);
        ti.flops = 0.0;
        ti.millis = 0.0;
        if (test_type == CSR) {
            csr_time(line, runs_num, &ti);
        } else if (test_type == HLL) {
            hll_time(line, runs_num, hack_size, &ti);
        }
        fprintf(off, "\"%s\",%f,%f\n", line, ti.millis, ti.flops);
    }
    fclose(iff);
    fclose(off);
}
