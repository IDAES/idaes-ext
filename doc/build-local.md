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
the scripts directory in order for patches (e.g. `ipopt.pc.patch`) to be
correctly applied.
- `ipopt.pc.patch` assumes that `coinmumps` and `coinhsl` are listed as
dependencies in the `ipopt.pc` file. If either is not found (i.e. we are
compiling without HSL), the patch application will fail.
- By default, Petsc is assumed to be installed in hard-coded, OS-dependent
location. If you are compiling locally, you have probably installed it
in your own preferred location. To override the default, set the `PETSC_DIR`
environment variable to the location where you have installed Petsc.
- If you have your own system-installed HSL libraries that are discoverable
by the linker (i.e. `LD_LIBRARY_PATH` is set), the COIN-OR solvers will likely
be able to find and link against them. In this case, you will see a
`HSL Present: NO` message (as there was no `coinhsl.zip` in the parent of the
working directory), but will likely still build HSL-enabled solvers.
