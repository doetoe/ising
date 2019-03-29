CC = g++ # clang++

all: ising

ising: ising.cpp
	$(CC) --std=c++14 -Wall -Wextra -Wpedantic -I. -o ising -O3 ising.cpp

ising.cpp: matrix.h

clean:
	rm -f *.o

realclean: clean
	rm -f ising

