#include "csr.h"
#include "mmio.h"
#include <string.h>

// Create the csr struct for generic data
int __csr_init(struct csr *csr, struct MatrixMarket *m) {
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

// Create the csr struct for int data
int csr_iinit(struct csr *csr, struct MatrixMarket *m)
{
    return __csr_init(csr, m);
}

// Create the csr struct for double data
int csr_dinit(struct csr *csr, struct MatrixMarket *m) {
    return __csr_init(csr, m);
}
