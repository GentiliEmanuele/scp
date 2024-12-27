#include "utils.h"
#include "mmio.h"
#include <stdio.h>
#include <stdlib.h>

static inline char *get_format_string(struct MatrixMarket *m) {
    if (mm_is_pattern(m->typecode)) {
        return "%d %d\n";
    } else {
        if (mm_is_integer(m->typecode)) {
            return "%d %d %ld\n";
        } else {
            return "%d %d %lf\n";
        }
    }
}

/**
 * read_mtx -- read a file in MatrixMarket format
 * @param path -- path of the file to read
 * @param mm   -- output struct
 * @return 1 if the file cannot be opened or the matrix is in a format not supported,
 * 0 otherwise.
 */
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

    // I  col indexes
    // J  row indexes
    // counter_nz real number of non zero items
    int counter_nz = 0;
    int *I = (int *) malloc(nz * sizeof(int));
    int *J = (int *) malloc(nz * sizeof(int));
    char *fmt = get_format_string(mm);
    int element_size = get_element_size(mm);
    void *val = malloc(nz * get_element_size(mm));
    for (int i = 0; i < nz; i++) {
        if (!mm_is_pattern(mm->typecode)) {
            if (mm_is_real(mm->typecode)) {
                fscanf(f, fmt, &I[i], &J[i], &((double*)val)[i]);
            } else {
                fscanf(f, fmt, &I[i], &J[i], &((int*)val)[i]);
            }
        } else {
            fscanf(f, fmt, &I[i], &J[i]);
        }
        if (mm_is_symmetric(mm->typecode) && (I[i] != J[i])) {
            counter_nz += 2;
        } else {
            counter_nz += 1;
        }
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }
    fclose(f);
    mm->num_rows = M;
    mm->num_cols = N;
    mm->nz = counter_nz;
    mm->data = val;
    mm->rows = J;
    mm->cols = I;
	return 0;
}