# Building IDAES Binaries

The build scripts are not intended for general use. IDAES users do not 
generally need to build IDAES binaries. 

Aside from macOS the Linux and Windows builds use Docker, and I'll assume 
we are using Docker Desktop.  The Linux ARM64 builds are done on macOS with 
Docker.  Since Windows Docker containers require Windows, the x86_64 builds
for Linux and Windows are done on Windows.

## Build Environments

On Windows, except for one excpetion noted below the git Bash shell can be used. 

### Docker

1. Install Docker Desktop.
2. Get the docker files from https://github.com/IDAES/idaes-ext/tree/main/docker
3. Go to the build-platform directory.
4. From the a bash shell run `sh build.sh {flavor}` Replace `{flavor}` with the 
  platform (el7, el8, ubuntu1804, ubuntu2004, ubuntu2204, or windows). On ARM64 
  there is no el7 or windows.
5. Usually testing is done by GitHub actions and the testing images are built on
  DockerHub, but ARM64 is not supported.  If you need to build testing containers,
  go to the test-platform directory and run `sh build.sh {platform}`.  The platform
  argument is slightly differnt than flavor.  For testing it indicates a specific
  Linux distribution (or Windows).

### macOS

1. Install the Xcode command line tools
2. Install cmake (cmake.org), and set it up to run from the command line
3. Install homebrew and a few things:
  * brew install gcc
  * brew install pkgconfig
  * brew install boost
  * brew install bash
  
You only need to rebuild PETSc when you want to update to a new version.  A
recommended location for this is $HOME/src.  For now, we compile without MPI, as
we currently aren't aiming to do super computing.

1. download https://petsc.org/release/download/
2. extract in ~/src
3. Go to PETSC source directory
4. Configure PETSc (-fPIC is so we can use the Metis and Mumps builds again)
   ./configure --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 \
      --download-metis --download-mumps --with-mumps-serial=1 \
      --prefix=$HOME/src/petsc-dist 
5. make
6. make install

## Building

### Linux

On Winodws make sure Docker Desktop is set to Linux mode.

1. Go to the build-extensions directory
2. Put coinhsl.zip in the extras directoy if you have it
3. In a bash shell (can use git bash on Windows) run 
  `sh build.sh {flavor} {arch}`, where arch is x86_64 or aarch64.

### Windows 

The `build.sh` script does not work for the Windows build. There is a PowerShell
script that will work. Make sure Docker Desktop is set to Windows mode (only works
on Windows).

1. Go to the build-extensions directory
2. Put coinhsl.zip in the extras directoy if you have it
3. In powershell run `.\build.ps1 windows --no-cache`

### MacOS

1. Checkout the idaes-ext repo.  
2. In the idaes-ext directory run the script to copy the
  source files to a new build location (
  `sh ./scripts/build_directory.sh ~/src/idaes-ext-build`).
3. Go to the build directory.
4. > bash ./scripts/compile_solvers.sh darwin
5. > bash ./scripts/compile_libs.sh darwin
6. > bash ./scripts/mac_collect.sh

## Testing

### x86 Windows and Linux

Testing is done with GitHub actions.  There are container images on docker hub.

### ARM64 Linux

There is a python script `docker_linux_tests.sh` in the idaes-ext/scripts directory. 
To run tests `sh docker_linux_tests.sh {platform}`. To help with debugging the
`test` container is left after the script completes.  You can run the container
interactively if needed.  Once done you will need to delete the container before
running another test.

### macOS (Apple Silicon Only)

Test on clean VM.

## Release Hashes

Collect all the tar files for a release in the same directory.

With IDAES installed, run ``idaes hash-release --path <path to tar files> --release <release>``
A text file with the hash of all the files will be written to the directory with the
release files.
