#include "utils.h"
#include "mmio.h"
#include <stdio.h>
#include <stdlib.h>

int read_mtx(const char *path, struct MatrixMarket *mm) {
    FILE* f = fopen(path, "r");
    if (f == NULL) {
        printf("Cannot open file %s\n", path);
        return 1;
    }
    
    if (mm_read_banner(f, &(mm->typecode)) != 0) {
        printf("Could not process Matrix Market banner.\n");
        return 1;
    }

    if (mm_is_complex(mm->typecode) && mm_is_matrix(mm->typecode) && mm_is_sparse(mm->typecode)) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(mm->typecode));
        return 1;
    }

    /* find out size of sparse matrix .... */
    int M, N, nz;
    int ret_code =  mm_read_mtx_crd_size(f, &M, &N, &nz);
    if (ret_code) {
        return 1;
    }

    int *I = (int *) malloc(nz * sizeof(int));
    int *J = (int *) malloc(nz * sizeof(int));
    double *val = (double *) malloc(nz * sizeof(double));
    for (int i = 0; i < nz; i++) {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }
    fclose(f);
    mm->num_rows = M;
    mm->num_cols = N;
    mm->nz = nz;
    mm->data = val;
    mm->rows = J;
    mm->cols = I;
	return 0;
}