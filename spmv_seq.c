#include "spmv_seq.h"

int d_spmv_csr_seq(double *res, struct csr *csr, double *v, int n) {
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns");
        return 1;
    }
    for (int i = 0; i < csr->num_rows; ++i) {
        res[i] = 0.0;
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            res[i] += ((double*)csr->data)[j] * v[csr->col_index[j]];
        }
    }
    return 0;
}

int i_spmv_csr_seq(int *res, struct csr *csr, int *v, int n) {
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns");
        return 1;
    }
    for (int i = 0; i < csr->num_rows; ++i) {
        res[i] = 0;
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            res[i] += ((int*)csr->data)[j] * v[csr->col_index[j]];
        }
    }
    return 0;
}