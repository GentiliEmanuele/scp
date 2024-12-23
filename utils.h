#include "mmio.h"

struct MatrixMarket {
    MM_typecode typecode;
    double *data;
    int    *rows;
    int    *cols;
    int    num_rows;
    int    num_cols;
    int    nz;
};

int read_mtx(const char *path, struct MatrixMarket *mm);