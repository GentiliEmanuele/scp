#ifndef CSR_H
#define CSR_H
#include "utils.h"
#include <stdlib.h>

typedef struct csr
{
    int *col_index;
    int *row_pointer;
    void *data;
    int num_rows;
    int num_cols;
} csr_t;

int csr_init(struct csr *csr, struct MatrixMarket *m);

void csr_cleanup(struct csr *csr);
#endif