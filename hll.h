#include "utils.h"

typedef struct hll
{
    int  *offsets;
    int  offsets_num;
    int  *col_index;
    void *data;
    int  data_num;
    int  hacks_num;
} hll_t;

int hll_init(struct hll *hll, int hack_size, struct MatrixMarket *mm);