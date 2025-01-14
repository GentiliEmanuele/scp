#include "spmv_openmp.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int d_spmv_csr_par(double *res, struct csr *csr, double *v, int n) {
    double *data = (double*)csr->data;
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns\n");
        return 1;
    }

    #pragma omp for schedule(static)
    for (int i = 0; i < csr->num_rows; ++i) {
        double sum = 0.0;
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            sum += data[j] * v[csr->col_index[j]];
        }
        res[i] = sum;
    }

    return 0;
}

int i_spmv_csr_par(int *res, struct csr *csr, int *v, int n) {
    int *data = (int*)csr->data;
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns\n");
        return 1;
    }

    #pragma omp for schedule(static)
    for (int i = 0; i < csr->num_rows; ++i) {
        int sum = 0;
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            sum += data[j] * v[csr->col_index[j]];
        }
        res[i] = sum;
    }
    return 0;
}

// Calculates number of rows of hack h
static inline int R(struct hll *hll, int h) {
    if (hll->hacks_num - 1 == h && hll->num_rows % hll->hack_size) {
        return hll->num_rows % hll->hack_size;
    }
    return hll->hack_size;
}

int d_spmv_hll_par(double *res, struct hll *hll, double *v, int n) {
    double *data = hll->data;
    if (n != hll->num_cols) {
        printf("matrix and vector must have the same number of columns\n");
        return 1;
    }

    int rows = R(hll, 0);
    #pragma omp for schedule(static)
    for (int h = 0; h < hll->hacks_num; ++h) {
        for (int r = 0; r < R(hll, h); ++r) {
            double sum = 0.0;
            for (int j = 0; j < hll->max_nzr[h]; ++j) {
                int k = hll->offsets[h] + r * hll->max_nzr[h] + j;
                sum += data[k] * v[hll->col_index[k]];
            }
            res[rows * h + r] = sum;
        }
    }
    return 0;
}

int i_spmv_hll_par(int *res, struct hll *hll, int *v, int n) {
    int *data = hll->data;
    if (n != hll->num_cols) {
        printf("matrix and vector must have the same number of columns\n");
        return 1;
    }

    int rows = R(hll, 0);
    #pragma omp for schedule(static)
    for (int h = 0; h < hll->hacks_num; ++h) {
        for (int r = 0; r < R(hll, h); ++r) {
            int sum = 0;
            for (int j = 0; j < hll->max_nzr[h]; ++j) {
                int k = hll->offsets[h] + r * hll->max_nzr[h] + j;
                sum += data[k] * v[hll->col_index[k]];
            }
            res[h * rows + r] = sum;
        }
    }
    return 0;
}