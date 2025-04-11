#include "spmv_seq.h"

int spmv_csr_seq(double *res, struct csr *csr, double *v, int n) {
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns");
        return 1;
    }
    for (int i = 0; i < csr->num_rows; ++i) {
        double sum = 0.0;
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            sum += csr->data[j] * v[csr->col_index[j]];
        }
        res[i] = sum;
    }
    return 0;
}

int spmv_hll_seq(double *res, struct hll *hll, double *v, int n) {
    if (n != hll->num_cols) {
        printf("matrix and vector must have the same number of columns\n");
        return 1;
    }
    int hack_size = hll->hack_size;
    for (int i = 0; i < hll->num_rows; ++i) {
        int hack = i / hack_size;
        int row_start = (i % hack_size) * hll->max_nzr[hack] + hll->offsets[hack];
        int row_end = row_start + hll->max_nzr[hack];
        double sum = 0.0;
        for (int j = row_start; j < row_end; ++j) {
            sum += hll->data[j] * v[hll->col_index[j]];
        }
        res[i] = sum;
    }
    return 0;  
}