#include "csr.h"
#include "spmv_openmp.h"
#include "utils.h"
#include "vec.h"
#include <stdio.h>
#include <cuda_runtime.h>

__global__ void product(double *res, int *row_pointer, double *data, int *col_index,  double *v, int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) {
        double sum = 0.0;
        for (int j = row_pointer[i]; j < row_pointer[i+1]; ++j) {
            sum += data[j] * v[col_index[j]];
        }
        res[i] = sum;
    }
}

void print_vec(double *v, int n) {
	for (int i = 0; i < n; i++) {
		printf("%d %lg\n", i, v[i]);
	}
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
    struct csr sm;
    if (csr_init(&sm, &mm)) {
        mtx_cleanup(&mm);
        return -1;
    }
    mtx_cleanup(&mm);
    int m = sm.num_rows;
    double *v = d_random(m);
    double *d_data;
    int *d_col_index;
    int *d_row_pointer;
    cudaError_t err = cuda_csr_init(&sm, &d_data, &d_col_index, &d_row_pointer);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        csr_cleanup(&sm);
        return -1;
    }
    double *d_result;
    err = cudaMalloc(&d_result, sm.num_rows * sizeof(double));
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_row_pointer);
    	csr_cleanup(&sm);
        return 1;
    }
    double *d_v;
    err = cudaMalloc(&d_v, sm.num_rows * sizeof(double));
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_row_pointer);
    	csr_cleanup(&sm);
    	cudaFree(d_result);
        return 1;        
    }
    err = cudaMemcpy(d_v, v, sm.num_rows * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_row_pointer);
    	csr_cleanup(&sm);
    	cudaFree(d_result);
        cudaFree(d_v);
        return 1;
    }
    // Perform SAXPY on 1M elements
    // 1 number of block in the grid
    // m number of the thread in the block
    product<<<2, 1024>>>(d_result, d_row_pointer, d_data, d_col_index, d_v, m);
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("error %d (%s): %s\n", err, cudaGetErrorName(err), cudaGetErrorString(err));
	    cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_row_pointer);
    	cudaFree(d_result);
        cudaFree(d_v);
    	csr_cleanup(&sm);
        return 1;
    }
    double *result = d_zeros(m);
    cudaMemcpy(result, d_result, sm.num_rows * sizeof(double), cudaMemcpyDeviceToHost);
    char r_path[256];
    sprintf(r_path, "%s.result", argv[1]);
    double *py_result = d_zeros(m);
    if (spmv_csr_par(py_result, &sm, v, m)) {
        printf("cannot execute csr product\n");
        cudaFree(d_data);
    	cudaFree(d_col_index);
    	cudaFree(d_row_pointer);
    	cudaFree(d_result);
    	csr_cleanup(&sm);
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
