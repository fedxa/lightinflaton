CC=gcc
CFLAGS=-std=c99
LOADLIBES=-lgsl -lgslcblas -lm

all: plotexample

plotexample.o chidecays.o: chidecays.h

plotexample: plotexample.o chidecays.o


clean:
	rm plotexample *.o
