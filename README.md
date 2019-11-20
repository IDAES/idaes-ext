# Building IDAES Binaries for Distribution

This is a script to compile solvers and function libraries for the IDAES biniary distribution.  
This script is just for IDAES team internal use.

## Setup a Build Environment

### Windows

Install MSYS2. MSYS2 provides a shell which will allow use of Linux style build tools. It also 
provides a convenient package manager (pacman) which allows for easy installation of build tools.

1. Go to https://www.msys2.org/
2. Download the x86_64 installer
3. Run the installer (the default options should be okay)

Open the MSYS2 MinGW 64-bit terminal (go to: start menu/MSYS2 64Bit/MSYS2 MinGW 64Bit).

Update the MSYS2 software:

```sh
pacman -Syu
```

Repeat the previous step until there are no more updates.

Install the build tools and libraries. Some packages installed are group packages, and 
pacman will prompt to select which packages you would like to install. Press "enter" for 
the default, which is all.:

```sh
pacman -S mingw-w64-x86_64-toolchain mingw-w64-x86_64-boost unzip patch make git zip
```

### Linux

We are currently building the Linux binaries with CentOS 7.  These generally seem to be
compatable with newer versions of Linux regardless of distribution.  We could build on
CentOS 6, but there are some minor compatablity issues with C++11, and you need to update
the git client to 2.x at least for the sript to work.

### OSX Build Environment

TBD

## Run the Build Script

```sh
git clone https://github.com/idaes/build-bin
```

If you have HSL library code copy the coinhsl directory the the ```build-bin``` directory.  
This will be copied to the Ipopt third party libraries before building.

```sh
cd build-bin
sh compile_libs.sh
sh compile_solvers.sh
```

Once the tar files are created, rename for the arciteture.  The files are formatted like ```idaes-*-{os}-{bits}.tar.gz```.  Where os is in {windows, linux, darwin} and bits is in {32, 64}.  For now
we are only building 64 bit.
