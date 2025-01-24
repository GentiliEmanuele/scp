#include "csr.h"
#include "hll.h"
#include "spmv_openmp.h"
#include "utils.h"
#include "vec.h"
#include <stdio.h>
#include <cuda_runtime.h>

#define num_of_rows(h, hack_size, hacks_num, num_rows) ((hacks_num - 1 == h && num_rows % hack_size) ? num_rows % hack_size : hack_size)

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

__global__ void cuda_spmv_hll(double *res, int hack_size, int hacks_num, double *data, int *offsets, int *col_index,  int *max_nzr, double *v, int n) {
    int h = blockDim.x * blockIdx.x + threadIdx.x;
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

void print_vec(double *v, int n) {
	for (int i = 0; i < n; i++) {
		printf("%d %lg\n", i, v[i]);
	}
}

cudaError_t cuda_hll_init(struct hll *hll, double *data, int *col_index, int *maxnzr, int *offsets) {
    cudaError_t err;
    err = cudaMalloc(&data, sizeof(double) * hll->data_num);
    if (err != cudaSuccess) {
        return err;
    }
    err = cudaMalloc(&col_index, sizeof(int) * hll->data_num);
    if (err != cudaSuccess) {
        cudaFree(data);
        return err;
    }
    err = cudaMalloc(&maxnzr, sizeof(int) * hll->hacks_num);
    if (err != cudaSuccess) {
        cudaFree(data);
        cudaFree(col_index);
        return err;
    }
    err = cudaMalloc(&offsets, sizeof(int) * hll->offsets_num);
    if (err != cudaSuccess) {
        cudaFree(data);
        cudaFree(col_index);
        cudaFree(maxnzr);
    }
    err = cudaMemcpy(data, hll->data, hll->data_num * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cudaFree(*data);
        cudaFree(*col_index);
        cudaFree(*row_pointer);
        return err;
    }
    err = cudaMemcpy(col_index, hll->col_index, hll->data_num * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cudaFree(data);
        cudaFree(col_index);
        cudaFree(row_pointer);
        return err;
    }
    err = cudaMemcpy(maxnzr, hll->max_nzr, hll->hacks_num * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cudaFree(data);
        cudaFree(col_index);
        cudaFree(row_pointer);
        return err;
    }
    err = cudaMemcpy(offsets, hll->offsets, hll->offsets_num * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cudaFree(data);
        cudaFree(col_index);
        cudaFree(row_pointer);
        cudaFree(maxnzr);
    }
    return err;
}

cudaError_t cuda_csr_init(struct csr *csr, double **data, int **col_index, int **row_pointer) {
    cudaError_t err;
    err = cudaMalloc(data, sizeof(double) * csr -> num_data);
    if (err != cudaSuccess) {
        return err;
    }
    err = cudaMalloc(col_index, sizeof(int) * csr ->num_data);
    if (err != cudaSuccess) {
        cudaFree(*data);
        return err;
    }
    err = cudaMalloc(row_pointer, sizeof(int) * (csr -> num_rows + 1));
    if (err != cudaSuccess) {
        cudaFree(*data);
        cudaFree(*col_index);
        return err;
    }
    err = cudaMemcpy(*data, csr -> data, csr->num_data * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cudaFree(*data);
        cudaFree(*col_index);
        cudaFree(*row_pointer);
        return err;
    }
    err = cudaMemcpy(*col_index, csr->col_index, csr->num_data * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cudaFree(*data);
        cudaFree(*col_index);
        cudaFree(*row_pointer);
        return err;
    }
    err = cudaMemcpy(*row_pointer, csr->row_pointer, (csr->num_rows + 1) * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cudaFree(*data);
        cudaFree(*col_index);
        cudaFree(*row_pointer);
        return err;
    }
    return cudaSuccess;
}

int main(int argc, char **argv) {
    if (--argc != 1) {
        printf("see usage: program matrix_path\n");
        return -1;
    }
    struct MatrixMarket mm;
    if (read_mtx(argv[1], &mm)) {
        return -1;
    }
    struct hll sm;
    if (hll_init(&sm, &mm)) {
        mtx_cleanup(&mm);
        return -1;
    }
    mtx_cleanup(&mm);
    int m = sm.num_rows;
    double *v = d_random(m);
    double *d_data;
    int *d_col_index;
    int *d_maxnzr;
    int *d_offsets;
    cudaError_t err = cuda_hll_init(&sm, &d_data, &d_col_index, &d_maxnzr, &d_offset);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        hll_cleanup(&sm);
        return -1;
    }
    double *d_result;
    err = cudaMalloc(&d_result, sm.num_rows * sizeof(double));
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_maxnzr);
    	hll_cleanup(&sm);
        return 1;
    }
    double *d_v;
    err = cudaMalloc(&d_v, sm.num_rows * sizeof(double));
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_maxnzr);
    	hll_cleanup(&sm);
    	cudaFree(d_result);
        return 1;        
    }
    err = cudaMemcpy(d_v, v, sm.num_rows * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_maxnzr);
    	hll_cleanup(&sm);
    	cudaFree(d_result);
        cudaFree(d_v);
        return 1;
    }
    // Perform SAXPY on 1M elements
    // 1 number of block in the grid
    // m number of the thread in the block
    cuda_spmv_hll<<<2, 1024>>>(d_result, sm.hack_size, sm.hacks_num, d_data, d_offsets, d_col_index, d_maxnzr, d_v, m);
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
	    cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_maxnzr);
    	cudaFree(d_result);
        cudaFree(d_v);
    	hll_cleanup(&sm);
        return 1;
    }
    double *result = d_zeros(m);
    cudaMemcpy(result, d_result, sm.num_rows * sizeof(double), cudaMemcpyDeviceToHost);
    double *py_result = d_zeros(m);
    if (spmv_hll_par(py_result, &sm, v, m)) {
        printf("cannot execute csr product\n");
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_maxnzr);
    	cudaFree(d_result);
    	hll_cleanup(&sm);
        return 1;
    }
    if (!d_veceq(result, py_result, m, 1e-6)) {
        printf("test failed\n");
    } else {
        printf("test passed\n");
    }
    printf("Result \n");
    print_vec(result, 10);
    printf("Pyresult \n");
    print_vec(py_result, 10);
    printf("v\n");
    print_vec(v, 10);
    free(py_result);
    free(result);
    free(v);
    cudaFree(d_data);
    cudaFree(d_col_index);
    cudaFree(d_row_pointer);
    cudaFree(d_result);
    csr_cleanup(&sm);
}
