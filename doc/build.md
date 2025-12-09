# Building IDAES Binaries

The build scripts are not intended for general use. IDAES users do not 
generally need to build IDAES binaries. 

Aside from macOS the Linux and Windows builds use Docker, and I'll assume 
we are using Docker Desktop.  The Linux ARM64 builds are done on macOS with 
Docker.  Since Windows Docker containers require Windows, the x86_64 builds
for Linux and Windows are done on Windows.

## Build Environments

On Windows, except for one exception noted below the git Bash shell can be used.

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
4. Manually set the HOMEBREW_PREFIX: `export HOMEBREW_PREFIX=$(/usr/local/bin/brew --prefix)`
5. Install required packages using the new `brew`: `/usr/local/bin/brew install gcc pkgconfig boost`
6. Download and install PETSc using the steps above.
   **NOTE**: You may need to manually set `CC`, `CCX`, `FC`, and `F77` and use the following:
   ```
   ./configure CC=$CC CXX=$CXX FC=$FC F77=$F77 --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 \
       --download-metis --download-mumps --with-mumps-serial=1 --prefix=~/repo/petsc-dist \
       --download-metis-cmake-arguments='-DCMAKE_POLICY_VERSION_MINIMUM=3.5' 
   ```

## Building

### Linux

TODO: Redo this whole section.

On Windows make sure Docker Desktop is set to Linux mode.

1. Go to the `build-scripts` directory
2. Put `coinhsl.tar.gz` in the `extras` directory if you have it
3. In a bash shell (can use git bash on Windows) run 
  `sh build.sh {flavor} {arch}`, where arch is x86_64 or aarch64.

### Windows

The `build.sh` script does not work for the Windows build. There is a PowerShell
script that will work. Make sure Docker Desktop is set to Windows mode (only works
on Windows).

1. Go to the `build-scripts` directory
2. Put `coinhsl.tar.gz` in the `extras` directory if you have it
3. In powershell run `.\build.ps1 windows --no-cache`

### MacOS

1. Checkout the idaes-ext repo.
2. If installing with coinhsl, move that tarball to peer with the idaes-ext repo.
3. In the idaes-ext directory run the script to copy the
  source files to a new build location (
  `sh ./scripts/build_directory.sh ~/repo/ext-build`).
4. Go to the build directory.
5. > bash ./scripts/compile_solvers.sh darwin # --with-hsl if using hsl
6. > bash ./scripts/compile_libs.sh darwin
7. > bash ./scripts/mac_collect.sh

## Release Hashes

Collect all the tar files for a release in the same directory.

With IDAES installed, run ``idaes hash-extensions --path <path to tar files> --release <release>``
A text file with the hash of all the files will be written to the directory with the
release files.

## Building without HSL

For testing purposes, it may be useful to build for multiple operating systems
(via Docker) without HSL. This may be accomplished by passing the `--without-hsl`
argument to `docker/build-extensions/build.sh`, after the OS name, architecture name,
repository URL, and repository branch.
This fifth positional argument to `build.sh` is interpreted as an argument (or arguments)
to send to `compile_solvers.sh`.
For example:
```bash
./build.sh ubuntu2204 x86_64 https://github.com/IDAES/idaes-ext.git main --without-hsl
```
