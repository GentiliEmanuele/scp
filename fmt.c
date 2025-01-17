#include "fmt.h"

void write_mrk_mtx(struct MatrixMarket *mm) {
    for (int i = 0; i < mm->nz; ++i) {
        printf("%d %d %lg", mm->rows[i], mm->cols[i], mm->data[i]);
    }
}

void write_csr_mtx(struct csr *csr, struct MatrixMarket *mm) {
    for (int i = 0; i < csr->num_rows; ++i) {
        int num_cols = csr->row_pointer[i+1] - csr->row_pointer[i];
        printf("row no. = %d has %d non-zero elements\n", i+1, num_cols);
    }
    for (int i = 0; i < csr->num_rows; ++i) {
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            int col = csr->col_index[j];
            double v = csr->data[j];
            printf("%d %d %lg\n", i, col, v);
        }
    }
}

void write_hll(struct hll *hll, struct MatrixMarket *mm) {
    printf("offsets:\n");
    for (int i = 0; i < hll->offsets_num; i++) {
        printf( i < hll->offsets_num - 1 ? "%d, " : "%d", hll->offsets[i]);
    }
    printf("\n");
    printf("indices:\n");
    int offp = 0;
    for (int i = 0; i < hll->data_num; ++i) {
        if (offp < hll->offsets_num && i == hll->offsets[offp]) {
            offp++;
            printf("|");
        }
        printf("%d ", hll->col_index[i]);
    }
    printf("|\n");
    printf("values:\n");
    offp = 0;
    for (int i = 0; i < hll->data_num; ++i) {
        if (offp < hll->offsets_num && i == hll->offsets[offp]) {
            offp++;
            printf("|");
        }
        printf("%lg", hll->data[i]);
    }
    printf("|\n");
}