#include "csr.h"
#include <string.h>

// Create the csr struct for generic data
int __csr_init(struct csr *csr, void *v, int *ii, int *jj, size_t b, int nz) {
    memcpy(csr->data, v, nz * b);
    memcpy(csr->col_index, jj, nz * b);
    csr->row_pointer[0] = 0;
    int prev = ii[0];
    int k = 1;
    for (int i = 0; i < nz; ++i) {
        int curr = ii[i];
        if (prev != curr) {
            prev = curr;
            csr->row_pointer[k++] = i;
        }
    }
    csr->row_pointer[k] = nz;
    return 0;
}

// Create the csr struct for int data
int csr_iinit(struct csr *csr, int *v, int *ii, int *jj, int nz)
{
    return __csr_init(csr, (void*)v, ii, jj, sizeof(int), nz);
}

// Create the csr struct for double data
int csr_dinit(struct csr *csr, double *v, int *ii, int *jj, int nz) {
    return __csr_init(csr, (void*)v, ii, jj, sizeof(double), nz);
}
