omp_test_csr:
	mkdir -p build
	gcc -g core/*.c omp/spmv_*.c omp/omp_test.c omp/omp_main_test_csr.c -o build/omp_test_csr.exe -fopenmp -lm -Icore -Iomp

omp_time_csr:
	mkdir -p build
	gcc -O3 core/*.c omp/spmv_*.c omp/omp_time.c omp/omp_main_time_csr.c -o build/omp_time_csr.exe -fopenmp -lm -Icore -Iomp

omp_test_hll:
	mkdir -p build
	gcc -g core/*.c omp/spmv_*.c omp/omp_test.c omp/omp_main_test_hll.c -o build/omp_test_hll.exe -fopenmp -lm -Icore -Iomp

omp_test_opt_hll:
	mkdir -p build
	gcc -g -Domp_hll_v2 core/*.c omp/spmv_*.c omp/omp_test.c omp/omp_main_test_hll.c -o build/omp_test_hll.exe -fopenmp -lm -Icore -Iomp

omp_time_hll:
	mkdir -p build
	gcc -O3 core/*.c omp/spmv_*.c omp/omp_time.c omp/omp_main_time_hll.c -o build/omp_time_hll.exe -fopenmp -lm -Icore -Iomp

omp_time_opt_hll:
	mkdir -p build
	gcc -O3 -Domp_hll_v2 core/*.c omp/spmv_*.c omp/omp_time.c omp/omp_main_time_hll.c -o build/omp_time_hll.exe -fopenmp -lm -Icore -Iomp

seq_time_csr:
	mkdir -p build
	gcc -O3 core/*.c omp/spmv_*.c omp/seq_time.c omp/seq_main_time_csr.c -o build/seq_time_csr.exe -fopenmp -lm -Icore -Iomp

seq_time_hll:
	mkdir -p build
	gcc -O3 core/*.c omp/spmv_*.c omp/seq_time.c omp/seq_main_time_hll.c -o build/seq_time_hll.exe -fopenmp -lm -Icore -Iomp

cuda:
	mkdir -p build/cuda
	cmake CMakeLists.txt -B build/cuda
	cd build/cuda
	make

footprint:
	mkdir -p build
	gcc -g info/mkinfo.c core/utils.c core/mmio.c core/hll.c core/csr.c -o build/mkinfo.exe -Icore

clean:
	rm -rf build/

