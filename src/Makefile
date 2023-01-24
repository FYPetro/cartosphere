CC = g++

INCLUDE = -I../include -I/usr/local/include/eigen3 -I/usr/local/include

CFLAGS = -g $(INCLUDE) -std=c++17

CPPFLAGS = -O3 -Wswitch

all: cartosphere

cartosphere: main.cpp functions.cpp mesh.cpp solver.cpp research.cpp globe.cpp
	pwd
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $+ -L/usr/local/lib -L/usr/local/lib -lfftw3 -largp

clean:
	-rm cartosphere

install:
	install -d $(prefix)/bin
	install -m 0755 cartosphere $(prefix)/bin

uninstall:
	-rm $(prefix)/bin/cartosphere

.PHONY: all clean check install uninstall