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
    wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.16.5.tar.gz; \
    tar -zxvf petsc-3.16.5.tar.gz; \
    mv petsc-3.16.5 petsc; \
    mkdir petsc-dist

RUN cd repo/petsc; \
    ./configure CPPFLAGS="-D_POSIX_C_SOURCE=199309L" --with-debug=0 --with-shared=0 --with-fortran-bindings=0 --download-mpich --download-mumps --download-hypre --download-sundials2 --download-superlu --download-scalapack --download-metis --download-parmetis --download-ptscotch --download-ml --download-suitesparse --download-strumpack --prefix=/repo/petsc-dist; \
    make;\
    make install;\
