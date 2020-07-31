FROM ubuntu:18.04

RUN apt-get -y update; \
    apt-get -y upgrade; \
    apt-get -y install ubuntu-dev-tools liblapack-dev libblas-dev unzip

RUN apt-get -y install gfortran libboost-all-dev cmake
RUN mkdir repo

ENTRYPOINT ["/bin/bash"]
