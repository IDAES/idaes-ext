FROM centos:centos8
MAINTAINER The CentOS Project <cloud-ops@centos.org>

RUN yum -y update; yum clean all
RUN yum -y install git make gcc gcc-c++ gcc-gfortran boost-devel cmake
RUN yum -y install wget unzip patch which file valgrind
RUN yum -y install flex bison
RUN yum -y install dnf-plugins-core
RUN yum config-manager --set-enabled powertools
RUN dnf --enablerepo=powertools -y install lapack-devel
RUN dnf -y install python3
RUN mkdir repo

RUN cd repo; \
    wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.15.0.tar.gz; \
    tar -zxvf petsc-3.15.0.tar.gz; \
    mv petsc-3.15.0 petsc; \
    mkdir petsc-dist; \
    cd petsc; \
    ./configure --with-debug=0 --with-shared=0 --with-fortran-bindings=0 --download-mpich --download-mumps --download-hypre --download-sundials --download-superlu --download-scalapack --download-metis --download-parmetis --download-ptscotch --download-ml --download-suitesparse --download-strumpack --prefix=/repo/petsc-dist; \
    make;\
    make install;\
