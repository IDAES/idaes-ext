FROM rockylinux:8

# Install dependencies
RUN yum -y update; yum clean all; \
    yum -y install git make cmake wget unzip patch which file; \
    yum -y install gcc gcc-c++ gcc-gfortran; \
    yum -y install dnf-plugins-core; \
    yum config-manager --set-enabled powertools; \
    yum -y install python3 lapack-devel findutils diffutils;

# Create directory structure
RUN mkdir repo

# Install PETSc
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
    wget https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.gz; \
    tar -zxvf boost_1_87_0.tar.gz; \
    cd ./boost_1_87_0; \
    sh bootstrap.sh --prefix=/usr/local/; \
    chmod a+x ./b2; \
    ./b2 install || true;

