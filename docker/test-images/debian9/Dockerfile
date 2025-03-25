FROM debian:9

RUN apt-get -y update; \
    ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get -y upgrade; \
    apt-get -y install libgfortran3 libgomp1 liblapack3 libblas3 git wget
RUN mkdir repo
RUN cd repo; \
    if [ `uname -m` = "aarch64" ]; then wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O ~/miniconda.sh; fi; \
    if [ `uname -m` = "x86_64" ]; then wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh; fi; \
    bash ~/miniconda.sh -b -p $HOME/miniconda
