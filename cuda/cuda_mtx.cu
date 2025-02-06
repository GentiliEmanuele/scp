#include "cuda_mtx.h"

cudaError_t cuda_csr_init(struct csr *csr, double **data, int **col_index, int **row_pointer) {
    cudaError_t err;
#ifdef CUDA_MANAGED
    err = cudaMallocManaged(data, sizeof(double) * csr->num_data);
#else
    err = cudaMalloc(data, sizeof(double) * csr->num_data);
#endif
    if (err != cudaSuccess) {
        printf("cannot allocate enough memory for data\n");
        return err;
    }
#ifdef CUDA_MANAGED
    err = cudaMallocManaged(col_index, sizeof(int) * csr->num_data);
#else
    err = cudaMalloc(col_index, sizeof(int) * csr->num_data);
#endif
    if (err != cudaSuccess) {
        printf("cannot allocate enough memory for col_index\n");
        cudaFree(*data);
        return err;
    }
#ifdef CUDA_MANAGED
    err = cudaMallocManaged(row_pointer, sizeof(int) * (csr->num_rows + 1));
#else
    err = cudaMalloc(row_pointer, sizeof(int) * (csr->num_rows + 1));
#endif
    if (err != cudaSuccess) {
        printf("cannot allocate enough memory for row_pointer\n");
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

cudaError_t cuda_hll_init(struct hll *hll, double **data, int **col_index, int **maxnzr, int **offsets) {
    cudaError_t err;
#ifdef CUDA_MANAGED
    err = cudaMallocManaged(data, sizeof(double) * hll->data_num);
#else
    err = cudaMalloc(data, sizeof(double) * hll->data_num);
#endif
    if (err != cudaSuccess) {
         printf("cannot allocate enough memory for data\n");
        return err;
    }
#ifdef CUDA_MANAGED
    err = cudaMallocManaged(col_index, sizeof(int) * hll->data_num);
#else
    err = cudaMalloc(col_index, sizeof(int) * hll->data_num);
#endif
    if (err != cudaSuccess) {
         printf("cannot allocate enough memory for col_index\n");
        cudaFree(*data);
        return err;
    }
#ifdef CUDA_MANAGED
    err = cudaMallocManaged(maxnzr, sizeof(int) * hll->hacks_num);
#else
    err = cudaMalloc(maxnzr, sizeof(int) * hll->hacks_num);
#endif
    if (err != cudaSuccess) {
        printf("cannot allocate enough memory for maxnzr\n");
        cudaFree(*data);
        cudaFree(*col_index);
        return err;
    }
#ifdef CUDA_MANAGED
    err = cudaMallocManaged(offsets, sizeof(int) * hll->offsets_num);
#else
    err = cudaMalloc(offsets, sizeof(int) * hll->offsets_num);
#endif
    if (err != cudaSuccess) {
        printf("cannot allocate enough memory for offsets\n");
        cudaFree(*data);
        cudaFree(*col_index);
        cudaFree(*maxnzr);
    }
    err = cudaMemcpy(*data, hll->data, hll->data_num * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cudaFree(*data);
        cudaFree(*col_index);
        return err;
    }
    err = cudaMemcpy(*col_index, hll->col_index, hll->data_num * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cudaFree(*data);
        cudaFree(*col_index);
        return err;
    }
    err = cudaMemcpy(*maxnzr, hll->max_nzr, hll->hacks_num * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cudaFree(*data);
        cudaFree(*col_index);
        return err;
    }
    err = cudaMemcpy(*offsets, hll->offsets, hll->offsets_num * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cudaFree(*data);
        cudaFree(*col_index);
        cudaFree(*maxnzr);
    }
    return err;
}