#ifndef SEQ_TIME_H
#define SEQ_TIME_H

typedef struct time_measurement_seq {
    double mean_time;
    double flops;
    double std_dev;
    double min;
    double max;
    int num_runs;
} time_measurement_seq_t;

int seq_time_csr(const char *file, int num_runs, time_measurement_seq_t *time_measurement);
int seq_time_hll(const char *file, int hack_size, int num_runs, time_measurement_seq_t *time_measurement);
void read_and_measure_seq_csr(char *path, int num_runs, char *out_path);
void read_and_measure_seq_hll(char *path, int hack_size, int num_runs, char *out_path);
#endif