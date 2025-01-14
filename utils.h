#ifndef UTILS_H
#define UTILS_H
#include "mmio.h"

struct MatrixMarket {
    MM_typecode typecode;
    void   *data;
    int    *rows; 
    int    *cols;
    int    num_rows;
    int    num_cols;
    int    nz;
};

static inline int get_element_size(struct MatrixMarket *m) {
    if (mm_is_integer(m->typecode)) {
        return sizeof(int);
    } else if (mm_is_real(m->typecode) || mm_is_pattern(m->typecode)) {
        return sizeof(double);
    } else {
        return 0;
    }
}

int read_mtx(const char *path, struct MatrixMarket *mm);
void read_and_measure_csr(char *path, int num_runs, int num_thread);
void read_and_measure_hll(char *path, int hack_size, int num_runs, int num_thread);
#endif