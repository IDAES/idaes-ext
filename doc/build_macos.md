The build scripts are not intended for general use. Most builds are done using
Docker containers. If you need to know how to set up a build environment for
Linux or Windows, see the Docker files.

To build the IDAES binary extensions on macOS follow the rough steps below.

Preliminary Setup
1. Install Xcode
2. Install the Xcode command line tools
3. Install cmake (cmake.org), and set it up to run from the command line
4. Install homebrew
  * brew install gcc
  * brew install pkgconfig
  * brew install boost
5. Build PETSc (no mpi)
  * download https://petsc.org/release/download/
  * extract in ~/src
  * ./configure --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 --download-metis --download-mumps --with-mumps-serial=1 --prefix=$HOME/src/petsc-dist
Build extensions
5. Get idaes-ext (clone or download from github)
6. If you have the HSL, put the coinhsl.zip file in the same directory as
  idaes-ext (not in it)
7. Go to the idaes-ext directory
8. Run build scripts.  The scripts assume the homebrew location is /opt/homebrew
  * > bash ./scripts/compile_solvers.sh
  * > bash ./scripts/compile_libs.sh

If the build completes without errors the *.tar.gz files will be in the
idaes-ext/dist-solvers, idaes-ext/dist-lib, and idaes-ext/dist-petsc directories.
To use these files with IDAES collect the tar file and extract them in the
directory given by the ``idaes bin-directory`` command. Once the build is
complete, it is recommended that you delete the idaes-ext directory.  
