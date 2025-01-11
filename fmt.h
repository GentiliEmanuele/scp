#include "csr.h"
#include "hll.h"
#include "utils.h"

void write_mrk_mtx(struct MatrixMarket *mm);

void write_csr_mtx(struct csr *csr, struct MatrixMarket *mm);

void write_hll(struct hll *hll, struct MatrixMarket *mm);