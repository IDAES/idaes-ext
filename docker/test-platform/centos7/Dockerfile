FROM centos:centos7
MAINTAINER The CentOS Project <cloud-ops@centos.org>

RUN yum -y update; yum clean all
RUN yum -y install lapack blas libgfortran libgomp git wget
RUN mkdir repo
RUN cd repo; \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh; \
    bash ~/miniconda.sh -b -p $HOME/miniconda
