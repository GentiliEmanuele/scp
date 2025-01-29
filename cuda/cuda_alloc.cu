#include <cuda_runtime.h>
#include <stdio.h>

int main(int argc, char *argv) {
    int n = 124406070/3;
    double *d1, *d2, *d3;
    cudaError_t err = cudaMalloc(&d1, n * sizeof(double));
    if (err != cudaSuccess) {
        printf("cannot allocate enough space for d1\n"):
    }
    cudaError_t err = cudaMalloc(&d2, n * sizeof(double));
    if (err != cudaSuccess) {
        printf("cannot allocate enough space for d2\n"):
    }
    cudaError_t err = cudaMalloc(&d3, n * sizeof(double));
    if (err != cudaSuccess) {
        printf("cannot allocate enough space for d3\n"):
    }
    cudaFree(*d1);
    cudaFree(*d2);
    cudaFree(*d3);
}