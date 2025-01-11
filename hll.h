#ifndef HLL_H
#define HLL_H
#include "utils.h"

typedef struct hll
{
    int  *offsets;
    int  offsets_num;
    int  *col_index;
    void *data;
    int  data_num;
    int  hacks_num;
    int  num_rows;
    int  num_cols;
    int *nzr;
    int *max_nzr;
} hll_t;

int hll_init(struct hll *hll, int hack_size, struct MatrixMarket *mm);
#endif