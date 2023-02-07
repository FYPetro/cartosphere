CXX = /usr/local/opt/llvm/bin/clang++

INCLUDE = -I$(ROOT_DIR)/include \
	-I/usr/local/include \
	-I/usr/local/opt/llvm/include \
	$(FFTWINC) $(EIGENINC)

CPPFLAGS = -g -std=c++17 -O3 -fopenmp $(INCLUDE)

LDFLAGS = -L/usr/local/opt/llvm/lib $(FFTWLIB)

all: cartosphere

cartosphere: main.cpp functions.cpp mesh.cpp solver.cpp research.cpp globe.cpp shapefile.cpp dsht.cpp
	pwd
	$(CXX) $(CPPFLAGS) -o $@ $+ $(LDFLAGS)

clean:
	-rm cartosphere

install:
	install -d $(prefix)/bin
	install -m 0755 cartosphere $(prefix)/bin

uninstall:
	-rm $(prefix)/bin/cartosphere

.PHONY: all clean check install uninstall
