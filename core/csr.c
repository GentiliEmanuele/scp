#include "csr.h"
#include "mmio.h"
#include <stdlib.h>
#include <string.h>

// Create the csr struct for generic data
int csr_init(struct csr *csr, struct MatrixMarket *m) {
    csr->num_rows = m->num_rows;
    csr->num_cols = m->num_cols;
    csr->num_data = m->nz;
    csr->data = malloc(m->nz * sizeof(double));
    if (csr->data == NULL) {
        return 1;
    }
    csr->row_pointer = malloc((csr->num_rows + 1) * sizeof(int));
    if (csr->row_pointer == NULL) {
        free(csr->data);
        return 1;
    }
    csr->col_index = malloc(m->nz * sizeof(int));
    if (csr->col_index == NULL) {
        free(csr->data);
        free(csr->row_pointer);
        return 1;
    }
    memcpy(csr->data, m->data, m->nz * sizeof(double));
    memcpy(csr->col_index, m->cols, m->nz * sizeof(int));
    csr->row_pointer[0] = 0;
    int prev = m->rows[0];
    int k = 1;
    for (int i = 0; i < m->nz; ++i) {
        int curr = m->rows[i];
        while (curr - prev > 1) {
            csr->row_pointer[k++] = i;
            ++prev;
        }
        if (curr - prev == 1) {
            prev = curr;
            csr->row_pointer[k++] = i;
        }
    }
    while (k <= m->num_rows) {
        csr->row_pointer[k++] = m->nz; 
    }
    return 0;
}

size_t csr_get_size(struct csr *csr, int nz) {
    return sizeof(double) * nz + sizeof(int) * nz + sizeof(int) * (csr->num_rows + 1) + sizeof(struct csr);
}

void csr_cleanup(struct csr *csr) {
    free(csr->col_index);
    free(csr->row_pointer);
    free(csr->data);
}