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

While MinGW does produce Windows native binaries, depending on linking options, some DLLs may 
be required. Add the MinWG/MSYS2 DLLs to your path. For example if MSYS2 was installed in the 
default location you would probably want to add C:\msys64\mingw64\bin.

### Linux

#### Ubuntu

On Windows, on the Ubuntu LTS quick install with Hyper-V. Otherwise Ubuntu LTS VM.

Setup the build environment:

```sh
sudo apt-get update
sudo apt-get install build-essential libboost-all-dev git gfortran
```

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
sh idaes_compile.sh
```


