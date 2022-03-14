#!/bin/sh

mkdir /c/repo
mkdir /c/repo/petc-dist
cd /c/repo
wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.16.5.tar.gz
tar -zxvf petsc-3.16.5.tar.gz
mv petsc-3.16.5 petsc
cd petsc
/usr/bin/python configure --with-debug=0 --with-shared=0 --with-mpi=1 --with-fortran-bindings=0 --download-mumps --download-scalapack --download-metis --download-suitesparse --prefix=/c/repo/petsc-dist
make
make install
