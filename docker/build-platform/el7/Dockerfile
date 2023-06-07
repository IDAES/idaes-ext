FROM centos:centos7

RUN yum -y update; yum clean all
RUN yum -y install git make gcc gcc-c++ gcc-gfortran valgrind
RUN yum -y install wget unzip patch which file lapack-devel blas-devel
RUN yum -y install flex bison
RUN yum -y install epel-release
RUN yum -y install cmake3
RUN ln -s /usr/bin/cmake3 /usr/bin/cmake
RUN yum -y remove git
RUN yum -y remove git-*
RUN yum -y install https://packages.endpointdev.com/rhel/7/os/x86_64/endpoint-repo.x86_64.rpm
RUN yum -y install git

RUN mkdir repo

RUN cd repo; \
    wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.19.2.tar.gz; \
    tar -zxvf petsc-3.19.2.tar.gz; \
    mv petsc-3.19.2 petsc; \
    mkdir petsc-dist

RUN cd repo; \
    wget https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz; \
    tar -zxvf boost_1_80_0.tar.gz; \
    cd ./boost_1_80_0; \
    sh bootstrap.sh --prefix=/usr/local/; \
    chmod a+x ./b2; \
    ./b2 install;

RUN cd repo/petsc; \
    ./configure CPPFLAGS="-D_POSIX_C_SOURCE=199309L" --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 --download-metis --download-mumps --with-mumps-serial=1 --prefix=/repo/petsc-dist; \
    make;\
    make install;
