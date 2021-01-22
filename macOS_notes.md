Compiling the IDAES binaries on macOS can be challenging. One of these biggest obstacles is that macOS ships with clang as the compiler, which includes a dated version of gcc. So far, we have had the best luck compiling these binaries with gcc version 9 or 10 available through homebrew (not clang/Xcode).

This document outlines the (approximate) steps to compile the solver and library binaries on a macOS machine. For more details, see https://github.com/IDAES/idaes-ext/issues/69.

# Step 0: Clone this repository onto your computer

# Step 1: Install homebrew and dependencies

The first step is to install homebrew, which is an unofficial package manager for macOS.  At the time of writing, do this by running the following command in the terminal:

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Check here for up to date instructions: https://brew.sh/

Next, we need to install missing packages which are dependencies for the script coinbrew (which greatly automates compiling coin-or solvers including Ipopt). Run this command:

```
brew install wget bash gcc dos2unix pkg-config cmake boost
```

See the instructions for more details: https://coin-or.github.io/coinbrew/

Tip: macOS ships with an old version of bash. You'll need the newer version to run coinbrew. It is easiest to install this with homebrew.

# Step 2: Determine the version of gcc

By default, homebrew will install gcc (and g++, gfortran) in /usr/local/bin. You'll likely find another version of gcc in /usr/bin that ships with macOS (or is installed with XCode). We now need to check which version of gcc we installed. In the terminal, run

```
ls -l /usr/local/bin/ | grep gcc
```

On my sytem, I see
```
c++-10 -> ../Cellar/gcc/10.2.0_2/bin/c++-10
cpp-10 -> ../Cellar/gcc/10.2.0_2/bin/cpp-10
g++-10 -> ../Cellar/gcc/10.2.0_2/bin/g++-10
gcc-10 -> ../Cellar/gcc/10.2.0_2/bin/gcc-10
```

This means I have gcc version 10 installed. We'll need this in a minute.

# Step 3: Locate where homebrew installed boost

The external functions need the C++ library boost. We need to find where homebrew installed the libraries. On my machine, they are in:

```
/usr/local/Cellar/boost/1.75.0_1/include
```

# Step 4: Edit the scripts and makefiles

In the file `/src/Makeline.in` do the following:
* Specify either `CC=gcc-9` or `CC=gcc-10`
* Specify either `CXX=g++-9` or `CXX=g++-10`
* Specify either `LINK=g++-9` or `LINK=g++-10`
* Specify `BOOST=` as the path to the boost header files

In the file `/scripts/compile_solvers.sh` do the following:
* Update the line `bash coinbrew build` line to specify the connect version of `gcc`, `g++`, and `gfortran`
* Update the line `cmake -DCMAKE_C_COMPILER=gcc-9 .` to `gcc-10` if needed.


In the file `/scripts/compile_libs.sh` do the following:
* Update the line `./configure` to specify the connect version of `gcc` and `gfortran`

# Step 5: Obtain coinhsl.zip and place one directory up from this directory

# Step 6: Compile

First run:

```
sh ./script/compile_solvers.sh darwin
```

Then run:

```
sh ./script/compile_libs.sh darwin
```

# Step 7: Copy and install files

Inside `/dist-lib`, you should find `idaes-lib-darwin-64.tar.gz`. Likewise, inside `/dist-solvers`, you should find `idaes-solvers-darwin-64.tar.gz`.

Create one folder called `binaries` somewhere on your computer and make note of the path. (When installing IDAES with the developer instructions, I put the folder `binaries` inside of the folder `ideas-pse`.)

Then to install the binaries, run:
```
idaes get-extensions --url file:./binaries/
```

You should replace `./binaries/` with the path to the folder containing the two `.tar.gz` files.

# Step 8: Run the IDAES test (for developer installation)

In the terminal, run:
```
pytest -m "not integration"
```

This will take ~5 minutes to run. The final output will include a summary of the tests (passed, failed, skipped, etc.)

# Tips
1. It is important to install all of the IDAES dependencies; this is handled when you pip install IDAES in a conda environment (see https://idaes-pse.readthedocs.io/en/stable/advanced_user_guide/advanced_install/index.html). Remember to switch to the conda environment where IDAES is installed to successfully compile the binaries.
