#include "spmv_seq.h"

int d_spmv_csr_seq(double *res, struct csr *csr, double *v, int n) {
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns");
        return 1;
    }
    for (int i = 0; i < csr->num_rows; ++i) {
        double sum = 0.0;
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            sum += ((double*)csr->data)[j] * v[csr->col_index[j]];
        }
        res[i] = sum;
    }
    return 0;
}

int i_spmv_csr_seq(int *res, struct csr *csr, int *v, int n) {
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns");
        return 1;
    }
    for (int i = 0; i < csr->num_rows; ++i) {
        int sum = 0;
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            sum += ((int*)csr->data)[j] * v[csr->col_index[j]];
        }
        res[i] = sum;
    }
    return 0;
}

int d_spmv_hll_seq(double *res, struct hll *hll, double *v, int n) {
    if (n != hll->num_cols) {
        printf("matrix and vector must have the same number of columns");
        return 1;
    }
    int z = 0;
    int start = 0;
    double *data = (double *)hll->data;
    double sum = 0.0;
    int k = 0;
    int i, j;
    for (i = 0; i < hll->num_rows; i++) {
        sum = 0.0;
        for (j = start; j - start < hll->nzr[i]; j++) {
            sum += data[j] * v[hll->col_index[j]];
        }
        res[i] = sum;
        start = hll->nzr[i];
    }
    return 0;
}