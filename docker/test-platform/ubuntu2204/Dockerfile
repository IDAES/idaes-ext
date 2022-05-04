FROM ubuntu:22.04

ENV DEBIAN_FRONTEND noninteractive

RUN mkdir repo
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime;
RUN apt-get -y update;
RUN apt-get -y upgrade;
RUN apt-get -y install libgfortran5 libgomp1 liblapack3 libblas3 git wget;
RUN cd repo; \
    if [ `uname -m` = "aarch64" ]; then wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O ~/miniconda.sh; fi; \
    if [ `uname -m` = "x86_64" ]; then wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh; fi; \
    bash ~/miniconda.sh -b -p $HOME/miniconda
