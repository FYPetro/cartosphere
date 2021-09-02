CC = g++

INCLUDE = -I../include -I../include/LAPACK -I/usr/local/include/eigen3 

CFLAGS = -g $(INCLUDE) -std=c++14

CPPFLAGS = -O3

all: cartosphere

cartosphere: main.cpp mesh.cpp
	pwd
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $+ -L/usr/local/lib -L/usr/local/lib/gcc/10 -lgfortran -lblas -llapack -llapacke -largp

clean:
	-rm cartosphere

install:
	install -d $(prefix)/bin
	install -m 0755 cartosphere $(prefix)/bin

uninstall:
	-rm $(prefix)/bin/cartosphere

.PHONY: all clean check install uninstall