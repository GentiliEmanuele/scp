#include "spmv_cuda.h"

#define num_of_rows(h, hack_size, hacks_num, num_rows) ((hacks_num - 1 == h && num_rows % hack_size) ? num_rows % hack_size : hack_size)

__global__ void cuda_spmv_hll(double *res, int hack_size, int hacks_num, double *data, int *offsets, int *col_index,  int *max_nzr, double *v, int n) {
    int h = blockDim.x * blockIdx.x + threadIdx.x;
    if (h < hacks_num) {
    	int rows = num_of_rows(0, hack_size, hacks_num, n);
    	for (int r = 0; r < num_of_rows(h, hack_size, hacks_num, n); ++r) {
        	double sum = 0.0;
        	for (int j = 0; j < max_nzr[h]; ++j) {
            		int k = offsets[h] + r * max_nzr[h] + j;
            		sum += data[k] * v[col_index[k]];
        	}
        	res[rows * h + r] = sum;
    	}
    }
}

__global__ void cuda_spmv_csr(double *res, int *row_pointer, double *data, int *col_index,  double *v, int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) {
        double sum = 0.0;
        for (int j = row_pointer[i]; j < row_pointer[i+1]; ++j) {
            sum += data[j] * v[col_index[j]];
        }
        res[i] = sum;
    }
}

__global__ void cuda_spmv_csr_v2(double *res, int *row_pointer, double *data, int *col_index,  double *v, int n) {
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    int warp_id = thread_id / 32;
    int lane = thread_id % 32;
    int row = warp_id;
    double sum = 0.0;
    if (row < n_rows) {
        int row_start = row_pointer[row];
        int row_end = row_pointer[row + 1];
        for (int element = row_start + lane; element < row_end; element += 32) {
            sum += data[element] * x[col_index[element]];
        }
    }
    sum = warp_reduce(sum);
    if (lane == 0 && row < n_rows) {
        res[row] = sum;
    }
}