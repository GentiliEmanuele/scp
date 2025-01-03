all:
	gcc -g csr.c mmio.c utils.c main.c hll.c spmv_seq.c -o main.exe

clean:
	rm main