CC = g++ # clang++

all: ising

ising: ising.cpp
	$(CC) --std=c++14 -Wall -I. -o ising -O3 ising.cpp

ising.cpp: matrix.h

isingtime: isingtime.cpp
	$(CC) --std=c++14 -I. -o isingtime -O3 isingtime.cpp

isingtime.cpp: matrix.h

clean:
	rm -f *.o

realclean: clean
	rm -f ising isingtime

