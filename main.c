#include "mmio.h"
#include "csr.h"
#include "hll.h"
#include "utils.h"
#include "spmv_seq.h"
#include "spmv_openmp.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void write_csr_mtx(struct csr *csr, struct MatrixMarket *mm) {
    for (int i = 0; i < csr->num_rows; ++i) {
        int num_cols = csr->row_pointer[i+1] - csr->row_pointer[i];
        printf("row no. = %d has %d non-zero elements\n", i+1, num_cols);
    }
    for (int i = 0; i < csr->num_rows; ++i) {
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            int col = csr->col_index[j]+1;
            if (mm_is_integer(mm->typecode)) {
                int v = ((int*)csr->data)[j];
                printf("%d %d %d\n", i+1, col, v);
            } else {
                double v = ((double*)csr->data)[j];
                printf("%d %d %f\n", i+1, col, v);
            }
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
        if (mm_is_integer(mm->typecode)) {
            printf("%d ", ((int*)hll->data)[i]);
        } else {
            printf("%f ", ((double*)hll->data)[i]);
        }
    }
    printf("|\n");
}

int main(int argc, char *argv[])
{
    struct MatrixMarket mm;
    if (read_mtx(argv[--argc], &mm)) {
        printf("Read error!");
        return 1;
    }
    struct hll sm;
    if (hll_init(&sm, 32, &mm)) { 
        printf("Error!\n");
        return 1;
    }

    printf("Matrix has %d rows and %d cols\n", sm.num_rows, sm.num_cols);

    srand(42);
    double *r = malloc(sm.num_cols * sizeof(double));
    double *v = malloc(sm.num_cols * sizeof(double));
/*
    clock_t start = clock();
    if (d_spmv_csr_par(r, &sm, v, sm.num_cols)) {
        return 1;
    }
    clock_t end = clock();
    printf("Total time spent (openmp): %f\n", ((double)(end - start) / CLOCKS_PER_SEC));
    start = clock();
    if (d_spmv_csr_seq(r, &sm, v, sm.num_cols)) {
        return 1;
    }
    end = clock();
    printf("Total time spent (serial): %f\n", ((double)(end - start) / CLOCKS_PER_SEC));
    */
    return 0;
}
