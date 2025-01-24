#ifndef HLL_H
#define HLL_H

#ifdef __cplusplus
extern "C"
{
#endif

#include "utils.h"
#include <stdlib.h>

typedef struct hll
{
    int  *offsets;
    int  offsets_num;
    int  *col_index;
    double *data;
    int  data_num;
    int  hacks_num;
    int  num_rows;
    int  num_cols;
    int  hack_size;
    int  *max_nzr;
} hll_t;

int hll_init(struct hll *hll, int hack_size, struct MatrixMarket *mm);

size_t hll_get_size(struct hll *hll);

void hll_cleanup(struct hll *hll);

#ifdef __cplusplus
} // extern "C"
#endif
#endif
