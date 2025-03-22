#ifndef SPMV_CUDA_H
#define SPMV_CUDA_H
#include <cuda_runtime.h>

#ifdef __cplusplus
extern "C"
{
#endif

__global__ void cuda_spmv_hll(double *res, int hack_size, int hacks_num, double *data, int *offsets, int *col_index,  int *max_nzr, double *v, int n);
__global__ void cuda_spmv_hll_v2(double *res, int hack_size, int hacks_num, double *data, int *offsets, int *col_index,  int *max_nzr, double *v, int n);
__global__ void cuda_spmv_hll_v3(double *res, int hack_size, int hacks_num, double *data, int *offsets, int *col_index,  int *max_nzr, double *v, int n);
__global__ void cuda_spmv_csr(double *res, int *row_pointer, double *data, int *col_index,  double *v, int n);
__global__ void cuda_spmv_csr_v2(double *res, int *row_pointer, double *data, int *col_index,  double *v, int n);
__global__ void cuda_spmv_csr_v3(double *res, int *row_pointer, double *data, int *col_index,  double *v, int n);
__global__ void cuda_spmv_csr_v4(double *res, int *row_pointer, double *data, int *col_index,  double *v, int n);

#ifdef __cplusplus
}
#endif
#endif
