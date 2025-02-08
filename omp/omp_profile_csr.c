#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "spmv_openmp.h"
#include "utils.h"
#include "vec.h"

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("usage:\n%s matrix threads_num\n", argv[0]);
        return 1;
    }
    struct MatrixMarket mm;
    if (read_mtx(argv[1], &mm)) {
        return 1;
    }
    long threads_num = strtol(argv[2], NULL, 10);
    struct csr sm;
    if (csr_init(&sm, &mm)) {
        mtx_cleanup(&mm);
        return 1;
    }
    
    int m = sm.num_rows;
    int n = sm.num_cols;
    double *r = d_zeros(m);
    srand(42);
    double *v = d_random(n);
    for (int i = 0; i < 50; i++) {
        if (spmv_csr_par(r, &sm, v, n)) {
            printf("warning: couldn't complete sparse matrix-vector product of run %d\n", i);
        }
    }
    csr_cleanup(&sm);
    mtx_cleanup(&mm);
    return 0;
}