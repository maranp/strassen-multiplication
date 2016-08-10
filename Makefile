CFLAGS=-Wextra -Wall -g -O3 -DNDEBUG -march=native
strass : strass.c main.c util.c

clean:
	rm -f strass *~
