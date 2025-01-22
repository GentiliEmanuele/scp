#include <stdio.h>
#include <cuda_runtime.h>

void gpuinfo(void) {
    int dev_count;
    cudaGetDeviceCount(&dev_count);
    for (int i = 0; i < dev_count; ++i) {
        cudaDeviceProp p;
        cudaGetDeviceProperties(&p, i);
        printf("#%d name=%s major=%d minor=%d\n", i, p.name, p.major, p.minor);
        printf("global memory=    %d\n", p.totalGlobalMem);
        printf("shared memory=    %d\n", p.sharedMemPerBlock);
        printf("constant memory=  %d\n", p.totalConstMem);
        printf("block registers=  %d\n", p.regsPerBlock);
        printf("warp size=        %d\n", p.warpSize);
        printf("threads per block=%d\n", p.maxThreadsPerBlock);
        printf("max block dim=    %d, %d, %d\n", p.maxThreadsDim[0], p.maxThreadsDim[1], p.maxThreadsDim[2]);
        printf("max grid dim=     %d, %d, %d\n", p.maxGridSize[0], p.maxGridSize[1], p.maxGridSize[2]);
    }
}

int main(int argc, char *argv[]) {
    gpuinfo();
}