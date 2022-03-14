FROM ubuntu:18.04

RUN apt-get -y update; \
    apt-get -y upgrade; \
    apt-get -y install ubuntu-dev-tools liblapack-dev libblas-dev unzip

RUN apt-get -y install gfortran libboost-all-dev cmake bison
RUN apt-get -y install flex valgrind
RUN mkdir repo

RUN cd repo; \
    wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.16.5.tar.gz; \
    tar -zxvf petsc-3.16.5.tar.gz; \
    mv petsc-3.16.5 petsc; \
    mkdir petsc-dist; \
    cd petsc; \
    ./configure --with-debug=0 --with-shared=0 --with-fortran-bindings=0 --download-mpich --download-mumps --download-hypre --download-sundials2 --download-superlu --download-scalapack --download-metis --download-parmetis --download-ptscotch --download-ml --download-suitesparse --download-strumpack --prefix=/repo/petsc-dist; \
    make;\
    make install;\
