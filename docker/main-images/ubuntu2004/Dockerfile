FROM ubuntu:20.04

# Install dependencies
RUN apt-get -y update; \
    ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get -y upgrade; \
    apt-get -y install ubuntu-dev-tools liblapack-dev libblas-dev unzip; \
    apt-get -y install gfortran cmake; \
    apt-get -y install pkg-config bison; \
    apt-get -y install flex valgrind;

# Create directory structure
RUN mkdir repo

# Install PETSc
# 3.20.6 is the newest version that works with the PETSc extras
# included in this repo
RUN cd repo; \
    wget https://gitlab.com/petsc/petsc/-/archive/v3.20.6/petsc-v3.20.6.tar.gz; \
    tar -zxvf petsc-v3.20.6.tar.gz; \
    mv petsc-v3.20.6 petsc; \
    mkdir petsc-dist; \
    cd petsc; \
    ./configure configure --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 --download-metis --download-mumps --with-mumps-serial=1 --prefix=/repo/petsc-dist; \
    make;\
    make install;

# Install boost
RUN cd repo; \
    wget https://boostorg.jfrog.io/artifactory/main/release/1.86.0/source/boost_1_86_0.tar.gz; \
    tar -zxvf boost_1_86_0.tar.gz; \
    cd ./boost_1_86_0; \
    sh bootstrap.sh --prefix=/usr/local/; \
    chmod a+x ./b2; \
    ./b2 install || true;
