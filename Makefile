all:
	gcc -g csr.c mmio.c utils.c main.c hll.c spmv_seq.c -o main.exe

openmp:
	gcc -g csr.c mmio.c omp_test.c omp_time.c utils.c hll.c spmv_openmp.c spmv_seq.c vec.c main.c -o main.exe -fopenmp -lm

clean:
	rm *.exe
