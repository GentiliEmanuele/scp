all:
	gcc -g csr.c mmio.c omp_test.c omp_time.c utils.c fmt.c hll.c spmv_openmp.c spmv_seq.c vec.c omp_main_test_csr.c -o main.exe -fopenmp -lm

openmp_test_csr:
	gcc -g csr.c mmio.c omp_test.c omp_time.c utils.c fmt.c hll.c spmv_openmp.c spmv_seq.c vec.c omp_main_test_csr.c -o openmp_test_csr.exe -fopenmp -lm

openmp_time_csr:
	gcc -g csr.c mmio.c omp_test.c omp_time.c utils.c hll.c fmt.c spmv_openmp.c spmv_seq.c vec.c omp_main_time_csr.c -o openmp_time_csr.exe -fopenmp -lm

openmp_test_hll:
	gcc -g csr.c mmio.c omp_test.c omp_time.c utils.c hll.c fmt.c spmv_openmp.c spmv_seq.c vec.c omp_main_test_hll.c -o openmp_test_hll.exe -fopenmp -lm

openmp_time_hll:
	gcc -g csr.c mmio.c omp_test.c omp_time.c utils.c hll.c fmt.c spmv_openmp.c spmv_seq.c vec.c omp_main_time_hll.c -o openmp_time_hll.exe -fopenmp -lm

openmp_pytest:
	gcc -g csr.c mmio.c utils.c hll.c omp_time.c spmv_seq.c spmv_openmp.c vec.c omp_main_pytest.c -o openmp_pytest.exe -fopenmp -lm

clean:
	rm *.exe
