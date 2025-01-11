all:
	gcc -g -Wall -Wextra -Wconversion csr.c mmio.c utils.c main.c hll.c spmv_seq.c -o main.exe

openmp:
	gcc -g -Wall -Wextra -Wconversion csr.c mmio.c utils.c hll.c spmv_openmp.c spmv_seq.c main.c -o main.exe -fopenmp -lm

clean:
	rm *.exe