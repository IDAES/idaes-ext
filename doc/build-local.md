# Building IDAES binaries locally

For testing purposes, it may be useful to build the binaries locally, i.e. not
using the Docker build scripts. This is also necessary for the MacOS builds
(see the MacOS section in `build.md` for more information).

You can build locally by running the following scripts:

- `scripts/compile_solvers.sh`: Compile solvers (including `k_aug`, `dot_sens`,
`petsc`, and all solver-associated libraries).
- `scripts/compile_libs.sh`: Compile shared libraries for external functions.

It is recommended to perform an "out-of-source" build. This can be done by
first running `scripts/build_directory.sh _build` from the top level of this
repository. This copies source code, makefiles, patch files, and build scripts
into the `_build` directory, where you will actually run the build.
Move into the build directory with `cd _build`.

The build scripts were not designed for local use, so there are a few gotchas:

- We look for patch files in a specific location relative to the original
working directory. I.e. you must run `compile_solvers.sh` from *one level above*
the scripts directory in order for patches (e.g. `CouenneMatrix.hpp.patch`) to be
correctly applied.
- By default, Petsc is assumed to be installed in hard-coded, OS-dependent
location. If you are compiling locally, you have probably installed it
in your own preferred location. To override the default, set the `PETSC_DIR`
environment variable to the location where you have installed Petsc.
- By default, we look for a `coinhsl.zip` file in the parent of the original
working directory. If this file is not found, we attempt to compile Ipopt with an
installed coinhsl library in a linker-accesible location. If no coinhsl library
can be found (or compiled from source in `coinhsl.zip`), the `compile_solvers.sh`
script will fail attempting to build `k_aug`, which cannot be built without HSL.
To avoid this failure, send the `--without-hsl` argument to `compile_solvers.sh`
after the OS name.
