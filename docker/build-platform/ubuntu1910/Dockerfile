FROM ubuntu:19.10

RUN apt-get -y update; \
    apt-get -y upgrade; \
    apt-get -y install ubuntu-dev-tools liblapack-dev libblas-dev unzip pkg-config

RUN apt-get -y install gfortran libboost-all-dev cmake
RUN mkdir repo

ENTRYPOINT ["/bin/bash"]
