all: main

main: main.c poly.c poly.h
	gcc -o main main.c poly.c poly.h -O3

clean:
	rm main *.o *.out *.gch