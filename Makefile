all: main

main: src/main.c src/poly.c src/poly.h
	gcc -o main src/main.c src/poly.c src/poly.h -O3

clean:
	rm main *.o *.out *.gch