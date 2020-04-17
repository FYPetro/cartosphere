CC = gcc

CFLAGS = -g -O0

all: cartosphere

cartosphere: main.c spherical.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $+

clean:
	-rm cartosphere

install:
	install -d $(prefix)/bin
	install -m 0755 cartosphere $(prefix)/bin

uninstall:
	-rm $(prefix)/bin/cartosphere

.PHONY: all clean check install uninstall