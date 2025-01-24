#include "csr.h"
#include "utils.h"
#include "vec.h"
#include <stdio.h>
#include <cuda_runtime.h>

__global__ void product(double *res, int *row_pointer, double *data, int *col_index,  double *v, int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) {
        double sum = 0.0;
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            sum += csr->data[j] * v[csr->col_index[j]];
        }
        res[i] = sum;
    }
}

void print_vec(double *v, int n) {
	for (int i = 0; i < n; i++) {
		printf("%d %lg\n", i, v[i]);
	}
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
    struct csr sm;
    if (csr_init(&sm, &mm)) {
        mtx_cleanup(&mm);
        return 1;
    }
    int m = sm.num_rows;
    int n = sm.num_cols;
    double *v = d_zeros(m);
    char v_path[256];
    sprintf(v_path, "%s.vector", argv[1]);
    read_vector(v, m, v_path);
    double *d_data;
    int *d_col_index;
    int *d_row_pointer;
    cudaMalloc(&d_data, sizeof(double) * sm.num_data);
    cudaMalloc(&d_col_index, sizeof(int) * sm.num_data);
    cudaMalloc(&d_row_pointer, sizeof(int) * (sm.num_rows + 1));
    double *d_result;
    cudaMalloc(&d_result, sm.num_rows * sizeof(double));
    double *d_v;
    cudaMalloc(&d_v, sm.num_rows * sizeof(double));
    cudaMemcpy(d_data, sm.data, sm.num_data * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_col_index, sm.col_index, sm.num_data * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_row_pointer, sm.row_pointer, (sm.num_rows + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, v, sm.num_rows * sizeof(double), cudaMemcpyHostToDevice);

    // Perform SAXPY on 1M elements
    // 1 number of block in the grid
    // m number of the thread in the block
    product<<<2, 1024>>>(d_result, d_row_pointer, d_data, d_col_index, d_v, m);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        goto cuda_error:
    }
    double *result = d_zeros(m);
    cudaMemcpy(result, d_result, sm.num_rows * sizeof(double), cudaMemcpyDeviceToHost);
    char r_path[256];
    sprintf(r_path, "%s.result", argv[1]);
    double *py_result = d_zeros(m);
    read_vector(py_result, m, r_path);
    if (!d_veceq(py_result, result, m, 1e-6)) {
        printf("test failed!\n");
    } else printf("test passed \n");
    print_vec(result, 10);
    print_vec(py_result, 10);
    free(py_result);
    free(result);
    free(v);
cuda_error:
    cudaFree(d_data);
    cudaFree(d_col_index);
    cudaFree(d_row_pointer);
    cudaFree(d_result);
    csr_cleanup(&sm);
    mtx_cleanup(&mm);
}
