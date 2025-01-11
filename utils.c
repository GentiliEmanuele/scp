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

static void readline(FILE *f, const char *fmt, int *row, int *col, void *val, int vp, struct MatrixMarket *mm) {
    if (!mm_is_pattern(mm->typecode)) {
        if (mm_is_real(mm->typecode)) {
            fscanf(f, fmt, col, row, &((double*)val)[vp]);            
        } else {
            fscanf(f, fmt, col, row, &((int*)val)[vp]);
        }
    } else {
        fscanf(f, fmt, col, row);
        ((double*)val)[vp] = 1.0;
    }
    // Adjust from one-based to zero-based
    *col -= 1;
    *row -= 1;
}

static int parse_rows_sy(FILE *f, struct MatrixMarket *mm) {
    struct vec3d *items = malloc(2 * mm->nz * sizeof(struct vec3d));
    void *data = malloc(mm->nz * get_element_size(mm));
    char *fmt = get_format_string(mm);
    int real_nz = 0;
    for (int k = 0; k < mm->nz; ++k) {
        readline(f, fmt, &items[real_nz].row, &items[real_nz].col, data, k, mm);
        items[real_nz].pos = k;
        int i = items[real_nz].row;
        int j = items[real_nz].col;
        ++real_nz;
        if (i != j) {
            items[real_nz].row = j;
            items[real_nz].col = i;
            items[real_nz].pos = k;
            ++real_nz;
        }
    }
    qsort(items, real_nz, sizeof(struct vec3d), qcmp);
    mm->rows = malloc(real_nz * sizeof(int));
    mm->cols = malloc(real_nz * sizeof(int));
    int element_size = get_element_size(mm);
    mm->data = malloc(real_nz * element_size);
    for (int i = 0; i < real_nz; ++i) {
        mm->rows[i] = items[i].row;
        mm->cols[i] = items[i].col;
        if (mm_is_real(mm->typecode) || mm_is_pattern(mm->typecode)) {
            memcpy(&((double*)mm->data)[i], &((double*)data)[items[i].pos], element_size);
        } else {
            memcpy(&((int*)mm->data)[i], &((int*)data)[items[i].pos], element_size);            
        }
    }
    free(data);
    printf("parsing symmetric matrix: number of non-zeros goes from %d to %d\n", mm->nz, real_nz);
    mm->nz = real_nz;
    return 0;
}


static int parse_rows_ns(FILE *f, struct MatrixMarket *mm) {
    int *rows = malloc(mm->nz * sizeof(int));
    int *cols = malloc(mm->nz * sizeof(int));
    void *values = malloc(mm->nz * get_element_size(mm));
    char *fmt = get_format_string(mm);
    for (int i = 0; i < mm->nz; i++) {
        readline(f, fmt, &rows[i], &cols[i], values, i, mm);
    }
    mm->rows = rows;
    mm->cols = cols;
    mm->data = values;
    return 0;
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

int *i_random(int n) {
    int *v = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        v[i] = rand();
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

int i_veceq(int *u, int *v, int n) {
    for (int i = 0; i < n; i++) {
        if (u[i] != v[i])
            return 1;
    }
    return 0;
}