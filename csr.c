#include "csr.h"
#include "mmio.h"
#include <stdlib.h>
#include <string.h>

// Create the csr struct for generic data
int csr_init(struct csr *csr, struct MatrixMarket *m) {
    csr->num_rows = m->num_rows;
    csr->num_cols = m->num_cols;
    csr->data = malloc(m->nz * get_element_size(m));
    if (csr->data == NULL) {
        return 1;
    }
    csr->row_pointer = malloc(csr->num_rows + 1);
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
    if (mm_is_pattern(m->typecode)) {
        double *data = (double*)m->data;
        for (int i = 0; i < m -> nz; i ++) {
            data[i] = 1;          
        }
    } else {
        memcpy(csr->data, m->data, m->nz * get_element_size(m));
    }
    memcpy(csr->col_index, m->cols, m->nz * sizeof(int));
    csr->row_pointer[0] = 0;
    int prev = m->rows[0];
    int k = 1;
    for (int i = 0; i < m->nz; ++i) {
        int curr = m->rows[i];
        if (prev != curr) {
            prev = curr;
            csr->row_pointer[k++] = i;
        }
    }
    csr->row_pointer[k] = m->nz;
    return 0;
}
