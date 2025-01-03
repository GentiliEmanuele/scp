#ifndef CSR_H
#define CSR_H
#include "utils.h"
#include <stddef.h>
typedef struct csr
{
    int *col_index;
    int *row_pointer;
    void *data;
    int num_rows;
    int num_cols;
} csr_t;

int csr_init(struct csr *csr, struct MatrixMarket *m);
#endif