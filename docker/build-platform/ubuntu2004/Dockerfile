FROM ubuntu:20.04

RUN apt-get -y update; \
    ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get -y upgrade; \
    apt-get -y install ubuntu-dev-tools liblapack-dev libblas-dev unzip

RUN apt-get -y install gfortran cmake
RUN apt-get -y install pkg-config bison
RUN apt-get -y install flex valgrind

RUN mkdir repo

RUN cd repo; \
    wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.18.1.tar.gz; \
    tar -zxvf petsc-3.18.1.tar.gz; \
    mv petsc-3.18.1 petsc; \
    mkdir petsc-dist; \
    cd petsc; \
    ./configure --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 --download-metis --download-mumps --with-mumps-serial=1 --prefix=/repo/petsc-dist; \
    make;\
    make install;

RUN cd repo; \
    wget https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz; \
    tar -zxvf boost_1_80_0.tar.gz; \
    cd ./boost_1_80_0; \
    sh bootstrap.sh --prefix=/usr/local/; \
    chmod a+x ./b2; \
    ./b2 install;
