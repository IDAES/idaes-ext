FROM centos:centos7

RUN yum -y update; yum clean all
RUN yum -y install git make gcc gcc-c++ gcc-gfortran boost-devel valgrind
RUN yum -y install wget unzip patch which file lapack-devel blas-devel
RUN yum -y install flex bison
RUN yum -y install epel-release
RUN yum -y install cmake3
RUN ln -s /usr/bin/cmake3 /usr/bin/cmake
RUN mkdir repo

RUN cd repo; \
    wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.17.1.tar.gz; \
    tar -zxvf petsc-3.17.1.tar.gz; \
    mv petsc-3.17.1 petsc; \
    mkdir petsc-dist

RUN cd repo/petsc; \
    ./configure CPPFLAGS="-D_POSIX_C_SOURCE=199309L" --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 --download-metis --download-mumps --with-mumps-serial=1 --prefix=/repo/petsc-dist; \
    make;\
    make install;\
    mkdir ../petsc-dist/include/mumps_libseq; \
    cp ./arch-*/externalpackages/git.mumps/libseq/mpi.h ../petsc-dist/include/mumps_libseq/;
