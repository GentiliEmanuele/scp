all:
	gcc -g csr.c mmio.c utils.c main.c hll.c -o main.exe

clean:
	rm main