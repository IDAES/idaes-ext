FROM ubuntu:22.04

RUN apt-get -y update; \
    ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get -y upgrade; \
    apt-get -y install libgfortran5 libgomp1 liblapack3 libblas3 git wget
RUN mkdir repo
RUN cd repo; \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh; \
    bash ~/miniconda.sh -b -p $HOME/miniconda
