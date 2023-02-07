CC = g++

INCLUDE = -I../include -I/usr/local/include -I/usr

CFLAGS = -g $(INCLUDE) -std=c++17

CPPFLAGS = -O3

all: cartosphere

cartosphere: main.cpp functions.cpp mesh.cpp solver.cpp research.cpp globe.cpp shapefile.cpp dsht.cpp
	pwd
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $+ \
		$(FFTWINC) $(FFTWLIB) $(EIGENINC)

clean:
	-rm cartosphere

install:
	install -d $(prefix)/bin
	install -m 0755 cartosphere $(prefix)/bin

uninstall:
	-rm $(prefix)/bin/cartosphere

.PHONY: all clean check install uninstall
