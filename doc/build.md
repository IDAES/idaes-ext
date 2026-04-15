# Building IDAES Binaries

The build scripts are **NOT** intended for general use. IDAES users do not
generally need to build IDAES binaries.

It is best to have both a Windows and a MacOS machine for this process.
The Linux builds (all platforms) can be completed on a MacOS, and
podman or Docker are both acceptable.
Windows + Docker is required for the Windows build.

## Build Environments

On Windows, except for one exception noted below, the git Bash shell can be used.

### CoinHSL

If you have the permissions to access the official HSL source, please make sure
that it is accessible in a `coinhsl.tar.gz` format. Place it in the
`docker/build-scripts/extras` directory.

### Docker

1. Install Docker Desktop.
3. Go to the `main-images` directory.
4. From the a bash shell run `sh build.sh {flavor}` Replace `{flavor}` with the 
  platform (el9, el8, ubuntu2404, ubuntu2004, ubuntu2204, or windows).
  - **NOTE**: If you are behind certain corporate firewalls/proxies, you will
    need to add lines to the beginning of the Dockerfiles for your SSL certificates.
  - **NOTE**: You can only build the Windows image if you are on a Windows machine.
    See [the official documentation](https://learn.microsoft.com/en-us/virtualization/windowscontainers/manage-docker/configure-docker-daemon).
  - **NOTE**: If you want to include `coinhsl.tar.gz`, place it in the `{flavor}`
    directory. It copies locally from there.

### macOS

There is logic in the `compile_solvers.sh` script to choose the correct
gcc, g++, etc., for Mac, based on `HOMEBREW_PREFIX`.

#### Apple Silicon

1. Install the Xcode command line tools
2. Install cmake (cmake.org), and set it up to run from the command line
3. Install homebrew and a few things:
  * brew install cmake gcc pkgconfig boost
  
You also need to install PETSc. The newest version that works with the custom
code in this repository is 3.20.6.

1. Download https://gitlab.com/petsc/petsc/-/archive/v3.20.6/petsc-v3.20.6.tar.gz
2. Extract somewhere appropriate for you (e.g., `mkdir ~/repo/` and extract there).
   **NOTE**: It is recommended rename the untarred directory to just `petsc`.
4. Go to PETSC source directory
5. Configure and build PETSc (-fPIC is so we can use the Metis and Mumps builds again)
   ```
   $ ./configure --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 \
       --download-metis --download-mumps --with-mumps-serial=1 --prefix=~/repo/petsc-dist \
       --download-metis-cmake-arguments='-DCMAKE_POLICY_VERSION_MINIMUM=3.5'
   $ make
   $ make install
   ```
   **NOTE**: If the configure is having difficulty downloading MUMPS for some reason,
   you can download it locally and point to it in the configure line instead using
   `--download-mumps=../MUMPS_5.6.1.tar.gz`
6. Set the `PETSC_DIR` environment var to wherever you installed PETSc.

#### Apple Intel

Apple has introduced a tool called [Rosetta](https://developer.apple.com/documentation/apple-silicon/about-the-rosetta-translation-environment)
to assist with the transition from Intel -> Silicon. We use this, while it's still available, for
building Apple Intel solvers.
Ensure to complete the steps below either in a different directory than from
Apple Silicon above **OR** delete old content from the Apple Silicon build.
The most common build issues occur when arm64 and x86_64 are mixed together
in your toolchain.

1. Install rosetta using `sudo softwareupdate --install-rosetta`
2. Create a new terminal and activate an x86_64 environment: `arch -x86_64 /bin/bash `
3. Install Homebrew for x86_64: `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"`.
   This will install brew in `/usr/local/bin/brew`.
4. Update the PATH: `export PATH=/usr/local/bin:$PATH`. This allows the correct `brew` to be found first.
5. Manually set the HOMEBREW_PREFIX: `export HOMEBREW_PREFIX=$(/usr/local/bin/brew --prefix)`
6. Install required packages using the new `brew`: `brew install gcc pkgconfig boost`
7. Download and install PETSc using the steps above.
   **NOTE**: You may need to manually set `CC`, `CCX`, `FC`, and `F77` and use the following:
   ```
   ./configure CC=$CC CXX=$CXX FC=$FC F77=$F77 --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 \
       --download-metis --download-mumps --with-mumps-serial=1 --prefix=~/repo/petsc-dist \
       --download-metis-cmake-arguments='-DCMAKE_POLICY_VERSION_MINIMUM=3.5' 
   ```
8. Set the `PETSC_DIR` environment var to wherever you installed PETSc.

## Release Hashes

Collect all the tar files for a release in the same directory.

Run `./scripts/hash-extensions release path-to-tar-files`.
A text file with the hash of all the tarballs will be written to the directory
in which the tarballs are located. They should then be copied over to the
`releases` directory.
