#include "utils.h"
#include "mmio.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <string.h>

struct vec3d {
    int row;
    int col;
    int pos;
};

static int inline lt(struct vec3d *a, struct vec3d *b) {
    if (a->row < b->row) {
        return 1;
    } else if ((a->row == b->row) && (a -> col < b -> col)) {
        return 1;
    }
    return 0;
}

static int qcmp(const void *a, const void *b) {
    return lt((struct vec3d*)a, (struct vec3d*)b);
}

static inline char *get_format_string(struct MatrixMarket *m) {
    if (mm_is_pattern(m->typecode)) {
        return "%d %d\n";
    } else {
        if (mm_is_integer(m->typecode)) {
            return "%d %d %ld\n";
        } else {
            return "%d %d %lg\n";
        }
    }
}

static void unpack(struct MatrixMarket *mm, int real_nz) {
    struct vec3d *items = malloc(real_nz * sizeof(struct vec3d));
    int p = 0;
    int *packed_rows = mm->rows;
    int *packed_cols = mm->cols;
    void *packed_values = mm->data;
    struct vec3d item;
    for (int k = 0; k < mm->nz; ++k) {
        int i = packed_rows[k];
        int j = packed_cols[k];
        item.row = i;
        item.col = j;
        item.pos = k;
        memcpy(&items[p], &item, sizeof(item));
        ++p;
        if (i != j) {
            item.row = j;
            item.col = i;
            item.pos = k;
            memcpy(&items[p], &item, sizeof(item));
            ++p;
        }
    }
    qsort(items, real_nz, sizeof(struct vec3d), qcmp);
    free(mm->rows);
    free(mm->cols);
    void *old_data = mm->data;
    mm->rows = malloc(real_nz * sizeof(int));
    mm->cols = malloc(real_nz * sizeof(int));
    int element_size = get_element_size(mm);
    mm->data = malloc(real_nz * element_size);
    for (int i = 0; i < real_nz; ++i) {
        mm->rows[i] = items[i].row;
        mm->cols[i] = items[i].col;
        if (mm_is_real(mm->typecode) || mm_is_pattern(mm->typecode)) {
            memcpy(&((double*)mm->data)[i], &((double*)mm->data)[items[i].pos], element_size);
        } else {
            memcpy(&((int*)mm->data)[i], &((int*)mm->data)[items[i].pos], element_size);            
        }
    }
    free(old_data);
    mm->nz = real_nz;
}

static int parse_rows(FILE *f, struct MatrixMarket *mm) {
    int *rows = malloc(mm->nz * sizeof(int));
    int *cols = malloc(mm->nz * sizeof(int));
    void *values = malloc(mm->nz * get_element_size(mm));
    char *fmt = get_format_string(mm);
    int real_nz = 0;
    for (int i = 0; i < mm->nz; i++) {
        if (!mm_is_pattern(mm->typecode)) {
            if (mm_is_real(mm->typecode)) {
                fscanf(f, fmt, &cols[i], &rows[i], &((double*)values)[i]);
            } else {
                fscanf(f, fmt, &cols[i], &rows[i], &((int*)values)[i]);
            }
        } else {
            fscanf(f, fmt, &cols[i], &rows[i]);
            ((double*)values)[i] = 1.0;
        }
        rows[i]--;  /* adjust from 1-based to 0-based */
        cols[i]--;
        real_nz++;
        if (mm_is_symmetric(mm->typecode) && rows[i] != cols[i]) {
            real_nz++;
        }
    }
    mm->rows = rows;
    mm->cols = cols;
    mm->data = values;
    if (mm->nz != real_nz && mm_is_symmetric(mm->typecode)) {
        printf("parsing a symmetric matrix, real number of non-zero values goes from %d to %d\n", mm->nz, real_nz);
        unpack(mm, real_nz);
    } else if (!mm_is_symmetric(mm->typecode) && mm->nz != real_nz) {
        printf("expected %d non-zeros values but got %d\n", mm->nz, real_nz);
        return 1;
    }
    return 0;
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
    return ir;
}

double *d_random(int n) {
    double *v = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        v[i] = (double)rand() / RAND_MAX;
    }
    return v;
    
}

int d_veceq(double *u, double *v, int n, double eps) {
    for (int i = 0; i < n; i++) {
        if (fabs(u[i] - v[i]) > eps)
            return 1;
    }
    return 0;
}