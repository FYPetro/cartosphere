CC = g++

incldir=../include

CFLAGS = -g -O0 -I$(incldir) -std=c++11

all: cartosphere

cartosphere: main.cpp mesh.cpp
	pwd
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $+ -L/usr/local/lib -largp

clean:
	-rm cartosphere

install:
	install -d $(prefix)/bin
	install -m 0755 cartosphere $(prefix)/bin

uninstall:
	-rm $(prefix)/bin/cartosphere

.PHONY: all clean check install uninstall