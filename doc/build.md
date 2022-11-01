# Building IDAES Binaries

The build scripts are not intended for general use. Most builds are done using
Docker containers.  IDAES users do not generally build IDAES binaries.

## Linux

### Building

Install docker desktop.

Copy the idaes-ext repo.

Build the docker build and test platforms.  The x86_64 Linux test images are also
built on DockerHub for use in GitHub actions.  To build the images just go to the
to either test-platform or build-platform in the docker directory then run the
build.sh (e.g. sh build.sh ubuntu1804) script with the flavor (OS of distro) to
build.  The available flavors are the subdirectories.  After running the build
script, a docker image should be available to build or test binaries.

If you want to compile with the HSL, put your coinhsl.zip file in the
./docker/build-extensions/extras directory.  If the HSL file is not available
the solvers will be built without it.

In the build-extensions directory run the build.sh script with the flavor to build
and the architecture {x86_64, aarch64} (e.g. sh build.sh ubuntu1804 aarch64). After
the build is complete tar files will be available in the ``tarballs`` subdirectory.

### Testing

Remove old builds from the test release on GitHub in the idaes-ext repo.  The test
containers will pull binaries from GitHub.

Trigger the GitHub actions tests or run the docker Linux test script with the
name of the platform to test.  

## Windows

### Building

On a Windows machine, install Docker Desktop.  Make sure docker is in Windows mode.
Use MSYS2 or the Git bash shell to run the build script in ./docker/build-platform
to create a windows build image (i.e. sh bash.sh windows).

Add coinhsl.zip to ./docker/build-extensions/extras if desired. Run the build.ps1
powershell script with the OS and arch names (i.e. build.ps1 windows x86_64).
The tarballs for the binaries will be in the tarballs directory.

### Testing

Remove old builds from the test release on GitHub in the idaes-ext repo.  The test
containers will pull binaries from GitHub.

Trigger the GitHub actions Windows test.

## macOS

To build the IDAES binary extensions on macOS follow the rough steps below.

### Preliminary Setup

1. Install Xcode
2. Install the Xcode command line tools
3. Install cmake (cmake.org), and set it up to run from the command line
4. Install homebrew
  * brew install gcc
  * brew install pkgconfig
  * brew install boost
  * brew install bash
  
### Build PETSc

You only need to rebuild PETSc when you want to update to a new version.  A
recommended location for this is $HOME/SRC.  For now, we compile without MPI, as
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

### Build

1. Get idaes-ext (clone or download from github)
2. If you have the HSL, put the coinhsl.zip file in the same directory as
  idaes-ext (not in it peer to).  Or in the directory were you plan to put the
  build directory later.
3. Go to the idaes-ext directory
  * (optional recommended) run the build directory script to create a build directory
    (e.g. sh ./scripts/build_directory.sh ~/src/idaes-ext-build)
  * Go to the build directory.
4. Run build scripts.  The scripts assume the homebrew location is /opt/homebrew, and
   homebrew is active.
  * > bash ./scripts/compile_solvers.sh darwin
  * > bash ./scripts/compile_libs.sh darwin
  * > bash ./scripts/mac_collect.sh

If the build completes without errors the *.tar.gz files will be in the
build directory.

### Test

Usually install the binaries on a minimal VM with IDAES-PSE installed and see
if the tests run properly.

## Release Hashes

Collect all the tar files for a release in the same directory.

With IDAES installed, run ``idaes hash-release --path <path to tar files> --release <release>``
A text file with the hash of all the files will be written to the directory with the
release files.
