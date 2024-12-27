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

int csr_iinit(struct csr *csr, struct MatrixMarket *m);
int csr_dinit(struct csr *csr, struct MatrixMarket *m);
#endif