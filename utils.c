#include "utils.h"
#include "mmio.h"
#include <stdio.h>
#include <stdlib.h>

int read_mtx(const char *path, double **v, int **I, int **J, int *M, int *N, int *nz)
{
    FILE* f = fopen(path, "r");
    if (f == NULL) {
        printf("Cannot open file %s\n", path);
        return 1;
    }
    
    MM_typecode matcode;
    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        return 1;
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        return 1;
    }

    /* find out size of sparse matrix .... */
    int ret_code =  mm_read_mtx_crd_size(f, M, N, nz);
    if (ret_code) {
        return 1;
    }

    int *ip = (int *) malloc(*nz * sizeof(int));
    int *jp = (int *) malloc(*nz * sizeof(int));
    double *val = (double *) malloc(*nz * sizeof(double));
    for (int i = 0; i < *nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &ip[i], &jp[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }
    *I = ip;
    *J = jp;
    *v = val;
    fclose(f);
      
	return 0;
}