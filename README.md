# Spherical Cartogram Library

Intended platforms:

- Windows 10/11 with Visual Studio (msvc19, msvc22)
- macOS (macOS 13 Ventura)

# Windows Instructions
## Installation
1. Install the newest Community edition of [Visual Studio](https://visualstudio.microsoft.com/ "Visual Studio")
2. Clone [vcpkg](https://vcpkg.io/en/docs/examples/installing-and-using-packages.html "vcpkg") into a desired location and bootstrap

	git clone https://github.com/microsoft/vcpkg.git
	.\vcpkg\bootstrap-vcpkg.bat

3. Add vcpkg folder to the `%PATH%` variable
4. Apply user-wide vcpkg root integration

	vcpkg integrate all

5. Clone this repository

	git clone https://github.com/FYPetro/cartosphere.git

## Build

1. Open .\msvc19\msvc19.sln
2. Build

# macOS Instructions

## System Configuration

1. Make sure command line tools are installed

	xcode-select --install

2. Install [Homebrew](https://brew.sh/ "Homebrew")
3. Install dependencies through brew

	brew install fftw[core,threads]
	brew install eigen
	brew install boost

4. Install C++ argparse

	git clone https://github.com/p-ranav/argparse
	cmake -DARGPARSE_BUILD_SAMPLES=on -DARGPARSE_BUILD_TESTS=on -DARGPARSE_INSTALL=on ..

5. You may need to add the following to the beginning of argparse's CMakeLists.txt

	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
	set(CMAKE_CXX_EXTENSIONS OFF)

6. Install llvm as our `-fopenmp`-supporting compiler

	brew install llvm

## Build

7. Clone this repository, compile, and install

	git clone https://github.com/FYPetro/cartosphere.git
	cd cartosphere
	make all
	make install
