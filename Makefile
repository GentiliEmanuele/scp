all:
	gcc -g csr.c mmio.c utils.c main.c -o main

clean:
	rm main