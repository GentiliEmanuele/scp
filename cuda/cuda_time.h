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

#ifdef __cplusplus
}
#endif
#endif