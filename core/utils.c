#include "utils.h"
#include "mmio.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#define IS_ZERO(x) (fabs(x) == 0.0)

struct vec3d {
    int row;
    int col;
    int pos;
};

static inline int lt(struct vec3d *a, struct vec3d *b) {
    if ((a->row < b->row) || (a->row == b-> row) && (a->col < b->col)){
        return -1;
    } else if ((a->row == b->row) && (a -> col == b -> col)) {
        return 0;
    } else {
        return 1;
    }
}

static inline int qcmp(const void *a, const void *b) {
    return lt((struct vec3d*)a, (struct vec3d*)b);
}

static int readline(FILE *f, struct vec3d *item, double *val, struct MatrixMarket *mm) {
    if (!mm_is_pattern(mm->typecode)) {
        fscanf(f, "%d %d %lg", &item->row, &item->col, val);
    } else {
        fscanf(f, "%d %d", &item->row, &item->col);
        *val = 1.0;
    }
    // Adjust from one-based to zero-based
    item->row -= 1;
    item->col -= 1;
    return !IS_ZERO(*val);
}

static int __parse_rows(FILE *f, struct MatrixMarket *mm) {
    int err;
    int is_symmetric = mm_is_symmetric(mm->typecode);
    struct vec3d *items = malloc((is_symmetric ? 2 : 1) * mm->nz * sizeof(struct vec3d));
    if (!items) {
        err = -1;
        goto no_items;
    }
    double *packed_data = malloc(mm->nz * sizeof(double));
    if (!packed_data) {
        err = -1;
        goto no_pck_data;
    }
    int real_nz = 0;
    int explicit_zeros = 0;
    for (int k = 0; k < mm->nz; k++) {
        if (readline(f, &items[real_nz], &packed_data[k], mm)) {
            int i = items[real_nz].row;
            int j = items[real_nz].col;
            items[real_nz].pos = k;
            ++real_nz;
            if (is_symmetric && i != j) {
                items[real_nz].row = j;
                items[real_nz].col = i;
                items[real_nz].pos = k;
                ++real_nz;
            }
        } else {
            ++explicit_zeros;
        }
    }
#ifdef SCP_VERBOSE
    if (is_symmetric)
        printf("symmetric matrix: number of non-zeros goes from %d to %d (explicit zeros = %d)\n", mm->nz, real_nz, explicit_zeros);
    else
        printf("non symmetric matrix: number of non-zeros %d (explicit zeros = %d)\n", mm->nz, explicit_zeros);
#endif
    qsort(items, real_nz, sizeof(struct vec3d), qcmp);
    int *rows = malloc(real_nz * sizeof(int));
    if (!rows) {
        err = -1;
        goto no_rows;
    }
    int *cols = malloc(real_nz * sizeof(int));
    if (!cols) {
        err = -1;
        goto no_cols;
    }
    double *data;
    if (mm_is_symmetric(mm->typecode)) {
        data = malloc(real_nz * sizeof(double));
        if (!data) {
            err = -1;
            goto no_data;
        }
        for (int i = 0; i < real_nz; ++i) {
            data[i] = packed_data[items[i].pos];
        }
    } else {
        data = packed_data;
    }
    for (int i = 0; i < real_nz; ++i) {
        rows[i] = items[i].row;
        cols[i] = items[i].col;
    }
    mm->rows = rows;
    mm->cols = cols;
    mm->data = data;
    mm->nz = real_nz;
    return 0;
no_data:
    free(cols);
no_cols:
    free(rows);
no_rows:
    free(packed_data);
no_pck_data:
    free(items);
no_items:
    return err;
}

static int parse_rows_sy(FILE *f, struct MatrixMarket *mm) {
    return __parse_rows(f, mm);
}


static int parse_rows_ns(FILE *f, struct MatrixMarket *mm) {
    return __parse_rows(f, mm);
}

static int parse_rows(FILE *f, struct MatrixMarket *mm) { 
    if (mm_is_symmetric(mm->typecode)) {
        return parse_rows_sy(f, mm);
    } else {
        return parse_rows_ns(f, mm);
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
        printf("cannot open file %s\n", path);
        return 1;
    }
    MM_typecode typecode;
    if (mm_read_banner(f, &(mm->typecode)) != 0) {
        printf("could not process Matrix Market banner\n");
        return 1;
    }

    if (mm_is_complex(mm->typecode) && mm_is_matrix(mm->typecode) && mm_is_sparse(mm->typecode)) {
        printf("sorry, this application does not support ");
        printf("matrix market type: [%s]\n", mm_typecode_to_str(mm->typecode));
        return 1;
    }

    int M, N, nz;
    int ret_code =  mm_read_mtx_crd_size(f, &M, &N, &nz);
    if (ret_code) {
        return 1;
    }

    mm->num_rows = M;
    mm->num_cols = N;
    mm->nz = nz;
    int ir = parse_rows(f, mm);
    fclose(f);
#ifdef SCP_VERBOSE
    int len = strlen(path);
    if (len < PATH_BUFFER_SZ) {
        strcpy(mm->__path_buffer, path);
        mm->__path_buffer[len] = 0;
        mm->path = mm->__path_buffer;
    } else {
        mm->path = malloc(len + 1);
        strcpy(mm->path, path);
    }
    printf("matrix %s:\n", path);
    printf("matrix has %d rows and %d cols and number of non-zeros %d\n", mm->num_rows, mm->num_cols, mm->nz);
#endif
    return ir;
}

void mtx_cleanup(struct MatrixMarket *mm) {
    free(mm->cols);
    free(mm->data);
    free(mm->rows);
#ifdef SCP_VERBOSE
    if (mm->path != mm->__path_buffer) {
        free(mm->path);
    }
#endif
}
