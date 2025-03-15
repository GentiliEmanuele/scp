#ifndef CUDA_TEST_H
#define CUDA_TEST_H
#ifdef __cplusplus
extern "C"
{
#endif

int csr_test(const char *path, int type);
int hll_test(char *path, int hack_size, int type);

#ifdef __cplusplus
}
#endif
#endif
