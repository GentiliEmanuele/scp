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
};

#define INIT_TIME_INFO(name) struct time_info name {.dev = 0, .millis = 0, .flops = 0};

float std_dev(float *samples, float avg, int n) {
    float dev = 0.0;
    for (int i = 0; i < runs_num; i++) {
        dev += (samples[i] - avg) * (samples[i] - avg);
    }
    return dev / n;
}
int csr_time(const char *path, int runs_num, struct time_info *ti);
int hll_time(const char *path, int runs_num, int hack_size, struct time_info *ti);

#ifdef __cplusplus
}
#endif
#endif