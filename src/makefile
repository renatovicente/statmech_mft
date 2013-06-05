CC = g++
OBJECTS = metropolis.o wanglandau.o

all: metropolis wanglandau

wanglandau : wanglandau.o
	$(CC)  wanglandau.o -lgsl -lgslcblas -lm -ligraph  -Wall -g -m64 -O2 -o wanglandau

metropolis : metropolis.o
	$(CC)  metropolis.o -lgsl -lgslcblas -lm -ligraph -Wall -g -m64 -O2 -o metropolis

wanglandau.o : wanglandau.cpp
metropolis.o : metropolis.cpp

clean :  $(OBJECTS)
	rm $(OBJECTS)
