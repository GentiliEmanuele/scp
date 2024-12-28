#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "csr.h"
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
    double *data = malloc(mm.nz * get_element_size(&mm));
    int *row_pointers = malloc((mm.num_rows + 1) * sizeof(int));
    int *col_index = malloc(mm.nz * sizeof(int));
    struct csr csr = {
        .col_index = col_index,
        .row_pointer = row_pointers,
        .data = data,
        .num_rows = mm.num_rows,
        .num_cols = mm.num_cols
    };
    if (csr_dinit(&csr, &mm)) { 
        printf("Error!\n");
    }
    write_csr_mtx(&csr, &mm);
    return 0;
}

