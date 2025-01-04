#include "spmv_openmp.h"
#include <omp.h>
#include <stdio.h>

int i_spmv_csr_par(int *res, struct csr *csr, int *v, int n) {
    int *data = (int*)csr->data;
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns");
        return 1;
    }
    int tid;
#pragma omp parallel private(tid)
{
    tid = omp_get_thread_num();
    printf("Thread id %d\n", tid);
    if (tid == 0) {
        printf("number of threads: %d\n", omp_get_num_threads());
    }
}
    int j;
    #pragma omp for private(j) schedule(static)
    for (int i = 0; i < csr->num_rows; ++i) {
        res[i] = 0;
        for (j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            res[i] += data[j] * v[csr->col_index[j]];
        }
    }
    return 0;
}

int d_spmv_csr_par(double *res, struct csr *csr, double *v, int n) {
    double *data = (double*)csr->data;
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns");
        return 1;
    }
    int tid;
#pragma omp parallel private(tid)
{
    tid = omp_get_thread_num();
    printf("Thread id %d\n", tid);
    if (tid == 0) {
        printf("number of threads: %d\n", omp_get_num_threads());
    }
}
    int j;
    #pragma omp for private(j) schedule(static)
    for (int i = 0; i < csr->num_rows; ++i) {
        res[i] = 0.0;
        for (j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            res[i] += data[j] * v[csr->col_index[j]];
        }
    }
    return 0;
}