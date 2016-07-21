CFLAGS=-Wall -g -O3 -DNDEBUG
strass : strass.c main.c util.c

clean:
	rm -f strass *~