#ifndef CSR_H
#define CSR_H
#include "utils.h"
#include <stdlib.h>

typedef struct csr
{
    int *col_index;
    int *row_pointer;
    double *data;
    int num_rows;
    int num_cols;
    int num_data;
} csr_t;

int csr_init(struct csr *csr, struct MatrixMarket *m);

size_t csr_get_size(struct csr *csr, int nz);

void csr_cleanup(struct csr *csr);
#endif