CFLAGS=-Wextra -Wall -g -O3 -DNDEBUG -march=native
LDFLAGS=-pthread
strass : strass.c main.c util.c

clean:
	rm -f strass *~
