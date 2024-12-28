#include "utils.h"

typedef struct hll
{
    void *hacks;
    int *col_index;
    int *offsets;
    void *data;
} hll_t;

int hll_iinit(struct hll *hll, struct MatrixMarket *mm, int hack_size);
int hll_dinit(struct hll *hll, struct MatrixMarket *mm, int hack_size);