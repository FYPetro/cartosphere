CXX = /usr/local/opt/llvm/bin/clang++

INCLUDE = -I$(ROOT_DIR)/include \
	-isystem/usr/local/include \
	-I/usr/local/opt/llvm/include \
	$(FFTWINC) $(EIGENINC)

SRC_FILES := $(wildcard *.cpp)

OBJ_FILES := $(patsubst %.cpp,%.o,$(SRC_FILES))

CPPFLAGS = -g -O2 -std=c++17 -fopenmp -D_RELEASE $(INCLUDE)

LDFLAGS = -L/usr/local/opt/llvm/lib \
	$(FFTWLIB) -lfftw3 -lfftw3_omp -lm \
	$(LOCALIB) -lglog

all: cartosphere

cartosphere: $(OBJ_FILES)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	pwd
	$(CXX) -c $(CPPFLAGS) -o $@ $+

clean:
	-rm cartosphere
	-rm *.o

install:
	install -d $(prefix)/bin
	install -m 0755 cartosphere $(prefix)/bin

uninstall:
	-rm $(prefix)/bin/cartosphere

.PHONY: all clean check install uninstall
