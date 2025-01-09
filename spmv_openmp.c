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
