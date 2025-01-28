#ifndef CUDA_MTX_H
#define CUDA_MTX_H
#include <cuda_runtime.h>
#include "csr.h"
#include "hll.h"

#ifdef __cplusplus
extern "C"
{
#endif

cudaError_t cuda_hll_init(struct hll *hll, double **data, int **col_index, int **maxnzr, int **offsets);
cudaError_t cuda_csr_init(struct csr *csr, double **data, int **col_index, int **row_pointer);

#ifdef __cplusplus
}
#endif
#endif