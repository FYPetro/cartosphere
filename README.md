# Cartosphere Library

License: `GPL-3.0-or-later`

`Cartosphere` is a `C++17` library that performs cartoglobe transformations natively on the unit sphere.
Based on [[Li2018](https://dx.doi.org/10.1080/15230406.2017.1408033)],
this current version `0.0.1` incorporates fixes documented in [[Li2023](https://www.proquest.com/openview/c8c54aced412bf581997e692f658d3df)].
The library currently only offers benchmarks,
and useful interfaces are planned.

Changelog: see `CHANGELOG.md`

Dependencies: `glog`, `fftw3`, `eigen3`, `omp`, `boost` (macOS)

Binaries may be compiled for:

- Windows 10/11 with Visual Studio (`msvc19`, `msvc22`)
- macOS (macOS 13 Ventura) with `llvm`

## Benchmark

```cartosphere benchmark```

## Instructions for Windows

1. Install the newest Community edition of [Visual Studio](https://visualstudio.microsoft.com/ "Visual Studio")
2. Clone [`vcpkg`](https://vcpkg.io/en/docs/examples/installing-and-using-packages.html "vcpkg") into a desired location and bootstrap
```
git clone https://github.com/microsoft/vcpkg.git
.\vcpkg\bootstrap-vcpkg.bat
```
3. Add `vcpkg` folder to the `%PATH%` variable
4. Apply user-wide `vcpkg` root integration
```
vcpkg integrate all
```
5. Clone this repository
```
git clone https://github.com/FYPetro/cartosphere.git
```
6. Open .\msvc19\msvc19.sln
7. Build

## Instructions for macOS

1. Install command-line tools after every major macOS upgrade:
```
xcode-select --install
```
2. Install [Homebrew](https://brew.sh/ "Homebrew — The Missing Package Manager for macOS (or Linux)")
3. Install dependencies through brew
```
brew install fftw[core,threads]
brew install eigen
brew install boost
brew install glog
```
4. Install C++ `argparse`
```
git clone https://github.com/p-ranav/argparse
cmake -DARGPARSE_BUILD_SAMPLES=on -DARGPARSE_BUILD_TESTS=on -DARGPARSE_INSTALL=on ..
```
5. You may need to add the following to the beginning of `argparse`'s `CMakeLists.txt`
```
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
```
6. Install `llvm` as our `-fopenmp`-supporting compiler
```
brew install llvm
```
7. Clone this repository, compile, and install
```
git clone https://github.com/FYPetro/cartosphere.git
cd cartosphere
make all
make install
```
