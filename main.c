#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "csr.h"
#include "hll.h"
#include "utils.h"

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
                printf("%d %d %f\n", i+1, col, v);
            } else {
                double v = ((double*)csr->data)[j];
                printf("%d %d %f\n", i+1, col, v);
            }
        }
    }
}

int main(int argc, char *argv[])
{
    struct MatrixMarket mm;
    if (read_mtx(argv[--argc], &mm)) {
        printf("Read error!");
        return 1;
    }
    struct hll hll;
    if (hll_dinit(&hll, &mm, 2)) { 
        printf("Error!\n");
    }
    return 0;
}

