omp_test_csr:
	mkdir -p build
	gcc -g core/*.c omp/spmv_*.c omp/omp_test.c omp/omp_main_test_csr.c -o build/omp_test_csr.exe -fopenmp -lm -Icore -Iomp

omp_time_csr:
	mkdir -p build
	gcc -g core/*.c omp/spmv_*.c omp/omp_time.c omp/omp_main_time_csr.c -o build/omp_time_csr.exe -fopenmp -lm -Icore -Iomp

omp_test_hll:
	mkdir -p build
	gcc -g core/*.c omp/spmv_*.c omp/omp_test.c omp/omp_main_test_csr.c -o build/omp_test_hll.exe -fopenmp -lm -Icore -Iomp

omp_time_hll:
	mkdir -p build
	gcc -g core/*.c omp/spmv_*.c omp/omp_time.c omp/omp_main_time_hll.c -o build/omp_time_csr.exe -fopenmp -lm -Icore -Iomp

omp_pytest:
	mkdir -p build
	gcc -g core/*.c omp/omp_main_pytest.c -o build/openmp_pytest.exe -fopenmp -lm -Icore -Iomp

info:
	mkdir -p build
	gcc -g core/mkinfo.c core/utils.c core/mmio.c core/hll.c core/csr.c -o build/mkinfo.exe -Icore

clean:
	rm build/*.exe
