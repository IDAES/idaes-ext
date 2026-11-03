FROM ubuntu:26.04

# Ensure commands fail on error
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Set noninteractive apt and timezone
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York

# Install dependencies
RUN apt-get update -y \
 && apt-get upgrade -y \
 && apt-get install -y \
      build-essential ca-certificates tzdata \
      ubuntu-dev-tools unzip cmake pkg-config wget \
      gfortran liblapack-dev libblas-dev libpython3-dev python3 python3-pip \
 && ln -fs /usr/share/zoneinfo/$TZ /etc/localtime \
 && dpkg-reconfigure -f noninteractive tzdata \
 && rm -rf /var/lib/apt/lists/*

# Create directory structure
RUN mkdir -p /repo

# Archives are extracted with Python's tarfile instead of GNU tar: the tar in
# Ubuntu 26.04 uses the openat2() syscall, which CPU emulators (QEMU/Rosetta,
# used for linux/amd64 builds on Apple Silicon) may not implement
RUN printf '%s\n' '#!/usr/bin/env python3' \
      'import sys, tarfile' \
      'with tarfile.open(sys.argv[1]) as t: t.extractall(filter="data")' \
      > /usr/local/bin/untar \
 && chmod +x /usr/local/bin/untar

# Copy the coinhsl.tar.gz file into the repo folder
# COPY ./hsl/coinhsl.tar.gz /repo/coinhsl.tar.gz

# Install PETSc
# Ubuntu 26.04 changes that affect PETSc 3.20.x:
#  * GCC 15 defaults to C23 (-std=gnu23). PETSc only gained GCC 15 fixes in
#    3.23.1, so pin C17 for PETSc
#  * Python 3.14 (two fixes, both needed by PETSc < 3.21):
#    - configure imports xdrlib, removed from the stdlib in 3.13
#    - config/gmakegen.py imports sysconfig._parse_makefile, which moved to
#      sysconfig.__main__ in 3.13
RUN python3 -m pip install --no-cache-dir --break-system-packages standard-xdrlib \
 && cd /repo \
 && wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.20.6.tar.gz \
 && untar petsc-3.20.6.tar.gz \
 && mv petsc-3.20.6 petsc \
 && mkdir petsc-dist \
 && cd petsc \
 && sed -i 's/^from sysconfig import _parse_makefile as parse_makefile$/try:\n  from sysconfig import _parse_makefile as parse_makefile\nexcept ImportError:  # Python >= 3.13\n  from sysconfig.__main__ import _parse_makefile as parse_makefile/' config/gmakegen.py \
 && export CMAKE_POLICY_VERSION_MINIMUM=3.5 \
 && ./configure --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 \
      --download-metis --download-mumps --with-mumps-serial=1 \
      CFLAGS="-std=gnu17" \
      --prefix=/repo/petsc-dist \
 && make \
 && make install

# Install boost
RUN cd /repo \
 && wget https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.gz \
 && untar boost_1_87_0.tar.gz \
 && cd boost_1_87_0 \
 && sh bootstrap.sh --prefix=/usr/local \
 && chmod a+x ./b2 \
 && ./b2 install
