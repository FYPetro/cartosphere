CC = g++

INCLUDE = -I../include 

CFLAGS = -g $(INCLUDE) -std=c++17 -Wswitch -Wclass-conversion

CPPFLAGS = -O3

all: cartosphere

cartosphere: main.cpp functions.cpp mesh.cpp solver.cpp research.cpp globe.cpp shapefile.cpp
	pwd
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $+ \
		$(FFTWINC) $(FFTWLIB) $(S2KITINC) $(S2KITOBJ) $(EIGENINC)

clean:
	-rm cartosphere

install:
	install -d $(prefix)/bin
	install -m 0755 cartosphere $(prefix)/bin

uninstall:
	-rm $(prefix)/bin/cartosphere

.PHONY: all clean check install uninstall
