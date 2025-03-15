#include "spmv_openmp.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int spmv_csr_par(double *res, struct csr *csr, double *v, int n, double *execution_time) {
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns\n");
        return 1;
    }
    double start = omp_get_wtime();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < csr->num_rows; ++i) {
        double sum = 0.0;
        #pragma omp simd reduction(+:sum)
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            sum += csr->data[j] * v[csr->col_index[j]];
        }
        res[i] = sum;
    }
    double end = omp_get_wtime();
    if (execution_time != NULL) {
        *execution_time = end - start;
    }
    return 0;
}

// Calculates number of rows of hack h
static inline int num_of_rows(struct hll *hll, int h) {
    if (hll->hacks_num - 1 == h && hll->num_rows % hll->hack_size) {
        return hll->num_rows % hll->hack_size;
    }
    return hll->hack_size;
}

int spmv_hll_par(double *res, struct hll *hll, double *v, int n, double *execution_time) {
    if (n != hll->num_cols) {
        printf("matrix and vector must have the same number of columns\n");
        return 1;
    }

    int rows = num_of_rows(hll, 0);
    double start = omp_get_wtime();
    #pragma omp parallel for schedule(static)
    for (int h = 0; h < hll->hacks_num; ++h) {
        for (int r = 0; r < num_of_rows(hll, h); ++r) {
            double sum = 0.0;
            #pragma omp simd reduction(+:sum)
            for (int j = 0; j < hll->max_nzr[h]; ++j) {
                int k = hll->offsets[h] + r * hll->max_nzr[h] + j;
                sum += hll->data[k] * v[hll->col_index[k]];
            }
            res[rows * h + r] = sum;
        }
    }
    double end = omp_get_wtime();
    if (execution_time != NULL) {
        *execution_time = end - start;
    }
    return 0;
}

int spmv_hll_par_v2(double *res, struct hll *hll, double *v, int n, double *execution_time) {
    if (n != hll->num_cols) {
        printf("matrix and vector must have the same number of columns\n");
        return 1;
    }
    int hack_size = hll->hack_size;
    double start = omp_get_wtime();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < hll->num_rows; ++i) {
        int hack = i / hack_size;
        int row_start = (i % hack_size) * hll->max_nzr[hack] + hll->offsets[hack];
        int row_end = row_start + hll->max_nzr[hack];
        double sum = 0.0;
        #pragma omp simd reduction(+:sum)
        for (int j = row_start; j < row_end; ++j) {
            sum += hll->data[j] * v[hll->col_index[j]];
        }
        res[i] = sum;
    }
    double end = omp_get_wtime();
    if (execution_time != NULL) {
        *execution_time = end - start;
    }
    return 0;  
}
