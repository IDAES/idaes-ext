FROM ubuntu:19.10

RUN apt-get -y update; \
    apt-get -y upgrade; \
    apt-get -y install libgfortran4 libgomp1 liblapack3 libblas3 git wget
RUN mkdir repo
RUN cd repo; \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh; \
    bash ~/miniconda.sh -b -p $HOME/miniconda

ENTRYPOINT ["/bin/bash"]
