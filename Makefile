all:
	gcc -g csr.c mmio.c utils.c main.c hll.c spmv_seq.c -o main.exe

openmp_test_csr:
	gcc -g csr.c mmio.c omp_test.c omp_time.c utils.c hll.c spmv_openmp.c spmv_seq.c vec.c omp_main_test_csr.c -o main.exe -fopenmp -lm

openmp_time_csr:
	gcc -g csr.c mmio.c omp_test.c omp_time.c utils.c hll.c spmv_openmp.c spmv_seq.c vec.c omp_main_time_csr.c -o main.exe -fopenmp -lm

openmp_test_hll:
	gcc -g csr.c mmio.c omp_test.c omp_time.c utils.c hll.c spmv_openmp.c spmv_seq.c vec.c omp_main_test_hll.c -o main.exe -fopenmp -lm

openmp_time_hll:
	gcc -g csr.c mmio.c omp_test.c omp_time.c utils.c hll.c spmv_openmp.c spmv_seq.c vec.c omp_main_time_hll.c -o main.exe -fopenmp -lm

clean:
	rm *.exe
