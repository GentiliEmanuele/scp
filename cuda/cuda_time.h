#ifndef CUDA_TIME_H
#define CUDA_TIME_H
#ifdef __cplusplus
extern "C"
{
#endif

struct time_info {
    float millis;
    float flops;
};

#define INIT_TIME_INFO(name) struct time_info name {.millis = 0, .flops = 0};

int csr_time(const char *path, int runs_num, struct time_info *ti);
int hll_time(const char *path, int runs_num, int hack_size, struct time_info *ti);

#ifdef __cplusplus
}
#endif
#endif