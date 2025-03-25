FROM rockylinux:9

# Install dependencies
RUN dnf -y update; dnf clean all; \
    dnf config-manager --set-enabled crb; \
    dnf -y install git make cmake wget unzip patch which file; \
    dnf -y install gcc gcc-c++ gcc-gfortran lapack-devel; \
    dnf -y install python3 python-devel findutils diffutils patchelf;

# Fix the symlinks: https://bugs.rockylinux.org/view.php?id=7855
RUN cd /usr/lib64; \
    ln -sf libblas.so.3.9.0 libblas.so.3; \
    ln -sf liblapack.so.3.9.0 liblapack.so.3;

# Create directory structure
RUN mkdir repo

# Copy the coinhsl.tar.gz file into the repo folder
# COPY ./hsl/coinhsl.tar.gz /repo/coinhsl.tar.gz

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


