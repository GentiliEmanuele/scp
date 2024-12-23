#include <stddef.h>
typedef struct csr
{
    int *col_index;
    int *row_pointer;
    void *data;
    int num_rows;
    int num_cols;
} csr_t;

int csr_iinit(struct csr *csr, int *v, int *ii, int *jj, int nz);
int csr_dinit(struct csr *csr, double *v, int *ii, int *jj, int nz);