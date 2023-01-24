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

3. Add vcpkg folder to the PATH variable
4. Install dependency and apply user-wide vcpkg root integration

	vcpkg install fftw3
	vcpkg install eigen3
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
3. Install dependencies

	brew install fftw
	brew install eigen

## Build
4. Clone this repository

	git clone https://github.com/FYPetro/cartosphere.git

5. Make

	cd cartosphere
	make all
	make install
