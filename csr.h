typedef struct csr
{
    int *col_index;
    int *row_pointer;
    double *data;
    int num_rows;
    int num_cols;
} csr_t;

double csr_get(csr_t *mtx, int i, int j);


