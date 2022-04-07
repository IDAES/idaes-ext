FROM ubuntu:22.04

RUN apt-get -y update; \
    export DEBIAN_FRONTEND=noninteractive; \
    ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime; \
    apt-get -y upgrade
Run apt-get -y install ubuntu-dev-tools liblapack-dev libblas-dev unzip
RUN apt-get -y install gfortran libboost-all-dev cmake
RUN apt-get -y install pkg-config bison
RUN apt-get -y install flex valgrind

RUN mkdir repo

RUN cd repo; \
    wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.16.5.tar.gz; \
    tar -zxvf petsc-3.16.5.tar.gz; \
    mv petsc-3.16.5 petsc; \
    mkdir petsc-dist; \
    cd petsc; \
    ./configure --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 --download-metis --download-mumps --with-mumps-serial=1 --prefix=/repo/petsc-dist; \
    make;\
    make install;\
