FROM ubuntu:20.04

# Ensure commands fail on error
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Set noninteractive apt and timezone
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York

# Install dependencies
RUN apt-get update -y \
 && apt-get install -y \
      ubuntu-dev-tools unzip cmake pkg-config wget \
      gfortran liblapack-dev libblas-dev libpython3-dev \
 && ln -fs /usr/share/zoneinfo/$TZ /etc/localtime \
 && dpkg-reconfigure -f noninteractive tzdata \
 && rm -rf /var/lib/apt/lists/*

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
