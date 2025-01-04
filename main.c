#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "csr.h"
#include "hll.h"
#include "utils.h"
#include "spmv_seq.h"
#include "spmv_openmp.h"

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
    struct csr sm;
    if (csr_init(&sm, &mm)) { 
        printf("Error!\n");
        return 1;
    }

    write_csr_mtx(&sm, &mm);

    int v[3] = {1,1,1};
    int r[3];
    if (i_spmv_csr_par(r, &sm, v, 3)) {
        return 1;
    }
    for (int i = 0; i < 3; i++) {
        printf("%d ", r[i]);
    }
    printf("\n");
    
    return 0;
}
