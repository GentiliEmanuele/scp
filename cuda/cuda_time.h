#ifndef CUDA_TIME_H
#define CUDA_TIME_H
#ifdef __cplusplus
extern "C"
{
#endif

struct time_info {
    float dev;
    float millis;
    float flops;
    float min;
    float max;
};

#define INIT_TIME_INFO(name) struct time_info name {.dev = 0, .millis = 0, .flops = 0, .min = 0, .max = 0};

int csr_time(const char *path, int runs_num, struct time_info *ti, int type);
int hll_time(const char *path, int runs_num, int hack_size, struct time_info *ti, int type);

#ifdef __cplusplus
}
#endif
#endif
