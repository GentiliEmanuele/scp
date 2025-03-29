#include "spmv_cuda.h"
#include <stdio.h>

#define num_of_rows(h, hack_size, hacks_num, num_rows) ((hacks_num - 1 == h && num_rows % hack_size) ? num_rows % hack_size : hack_size)


__global__ void cuda_spmv_hll_v0(double *res, int hack_size, int hacks_num, double *data, int *offsets, int *col_index,  int *max_nzr, double *v, int n) {
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

__global__ void cuda_spmv_hll_v1(double *res, int hack_size, int hacks_num, double *data, int *offsets, int *col_index,  int *max_nzr, double *v, int n) {
    int i = blockDim.x *blockIdx.x + threadIdx.x;
    if (i < n) {
        int hack = i / hack_size;
        int row_start = (i % hack_size) * max_nzr[hack] + offsets[hack];
        int row_end = row_start + max_nzr[hack];
        double sum = 0.0;
        for (int j = row_start; j < row_end; ++j) {
            sum += data[j] * v[col_index[j]];
        }
        res[i] = sum;
    }
}

__global__ void cuda_spmv_hll_v2(double *res, int hack_size, int hacks_num, double *data, int *offsets, int *col_index,  int *max_nzr, double *v, int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    extern __shared__ double vTile[];
    if (threadIdx.x < n) vTile[threadIdx.x] = v[threadIdx.x];
    if (i < n) {
        int hack = i / hack_size;
        int row_start = (i % hack_size) * max_nzr[hack] + offsets[hack];
        int row_end = row_start + max_nzr[hack];
        double sum = 0.0;
        for (int j = row_start; j < row_end; ++j) {
	    if (col_index[j] < 1024) {
		sum += data[j] * vTile[col_index[j]];
	    } else {
	        sum += data[j] * v[col_index[j]];
	    }
        }
	res[i] = sum;
    }

}
__global__ void cuda_spmv_hll_v3(double *res, int hack_size, int hacks_num, double *data, int *offsets, int *col_index,  int *max_nzr, double *v, int n) {
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    int warp_id = thread_id / 32;
    int lane = thread_id % 32;
    int row = warp_id;
    double sum = 0.0;
    if (row < n) {
        int hack = row / hack_size;
        int row_start = (row % hack_size) * max_nzr[hack] + offsets[hack];
        int row_end = row_start + max_nzr[hack];
        for (int element = row_start + lane; element < row_end; element += 32) {
            sum += data[element] * v[col_index[element]];
        } 
        for (int offset = 16; offset > 0; offset /= 2) {
            sum += __shfl_down_sync(0xFFFFFFFF, sum, offset);
        }
        if (lane == 0) res[row] = sum;
    }
}

__global__ void cuda_spmv_hll_v4(double *res, int hack_size, int hacks_num, double *data, int *offsets, int *col_index,  int *max_nzr, double *v, int n) {
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    int warp_id = thread_id / 32;
    int lane = thread_id % 32;
    int row = warp_id;
    double sum = 0.0;
    extern __shared__ double results[];
    if (row < n) {
        int hack = row / hack_size;
        int row_start = (row % hack_size) * max_nzr[hack] + offsets[hack];
        int row_end = row_start + max_nzr[hack];
	results[threadIdx.x] = 0.0;
        for (int element = row_start + lane; element < row_end; element += 32) {
            results[threadIdx.x] += data[element] * v[col_index[element]];
        } 
	sum = results[threadIdx.x];
	for (int offset = 16; offset > 0; offset /= 2) {
            sum += __shfl_down_sync(0xFFFFFFFF, sum, offset);
        }
        if (lane == 0) res[row] = sum;
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
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    extern __shared__ double vTile[];
    if (threadIdx.x < n) vTile[threadIdx.x] = v[threadIdx.x];
    if (i < n) {
	double sum = 0.0;
        for (int j = row_pointer[i]; j < row_pointer[i+1]; ++j) {
	    if (col_index[j] < 1024) {
		sum += data[j] * vTile[col_index[j]];
	    } else {
		sum += data[j] * v[col_index[j]];
	    }
        }
	res[i] = sum;
    }
}


__global__ void cuda_spmv_csr_v3(double *res, int *row_pointer, double *data, int *col_index,  double *v, int n) {
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    int warp_id = thread_id / 32;
    int lane = thread_id % 32;
    double sum = 0.0;
    if (warp_id > n) return;
    int row_start = row_pointer[warp_id];
    int row_end = row_pointer[warp_id + 1];
    for (int element = row_start + lane; element < row_end; element += 32) {
	sum += data[element] * v[col_index[element]];
    }
    for (int offset = 16; offset > 0; offset /= 2) {
        sum += __shfl_down_sync(0xFFFFFFFF, sum, offset);
    }
    if (lane == 0) res[warp_id] = sum;
}


__global__ void cuda_spmv_csr_v4(double *res, int *row_pointer, double *data, int *col_index,  double *v, int n) {
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    int warp_id = thread_id / 32;
    int lane = thread_id % 32;
    double sum = 0.0;
    extern __shared__ double results[];
    if (warp_id < n) {
	int row_start = row_pointer[warp_id];
    	int row_end = row_pointer[warp_id + 1];
        results[threadIdx.x] = 0.0;
    	for (int element = row_start + lane; element < row_end; element += 32) {
		results[threadIdx.x] += data[element] * v[col_index[element]];
    	}
	sum = results[threadIdx.x];
    	for (int offset = 16; offset > 0; offset /= 2) {
        	sum += __shfl_down_sync(0xFFFFFFFF, sum, offset);
    	}
    	if (lane == 0) res[warp_id] = sum;
    }
}

