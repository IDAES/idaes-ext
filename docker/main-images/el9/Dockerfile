FROM rockylinux:9

# Ensure commands fail on error
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Install dependencies
RUN dnf -y update \
 && dnf -y install dnf-plugins-core \
 && dnf config-manager --set-enabled crb \
 && dnf -y install epel-release \
 && dnf -y install \
        git make cmake wget unzip patch which file \
        gcc gcc-c++ gcc-gfortran lapack-devel \
        python3 python3-devel findutils diffutils patchelf \
 && dnf clean all

# Fix LAPACK/BLAS symlinks
# See: https://bugs.rockylinux.org/view.php?id=7855
RUN cd /usr/lib64 \
 && ln -sf libblas.so.3.9.0 libblas.so.3 \
 && ln -sf liblapack.so.3.9.0 liblapack.so.3

# Create directory structure
RUN mkdir -p /repo

# Copy the coinhsl.tar.gz file into the repo folder
# COPY ./coinhsl.tar.gz /repo/coinhsl.tar.gz

# Install PETSc
RUN cd /repo \
 && wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.20.6.tar.gz \
 && tar -zxvf petsc-3.20.6.tar.gz \
 && mv petsc-3.20.6 petsc \
 && mkdir petsc-dist \
 && cd petsc \
 && ./configure --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 \
      --download-metis --download-mumps --with-mumps-serial=1 \
      --prefix=/repo/petsc-dist \
 && make \
 && make install

# Install boost
RUN cd /repo \
 && wget https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.gz \
 && tar -zxvf boost_1_87_0.tar.gz \
 && cd boost_1_87_0 \
 && sh bootstrap.sh --prefix=/usr/local \
 && chmod a+x ./b2 \
 && ./b2 install
