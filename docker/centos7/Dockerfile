FROM centos:centos7
MAINTAINER The CentOS Project <cloud-ops@centos.org>

RUN yum -y update; yum clean all
RUN yum -y install git make gcc gcc-c++ gcc-gfortran boost-devel
RUN yum -y install wget

##
# Pull down the repo
##
RUN mkdir repo; \
    cd repo; \
    git clone https://github.com/dangunter/idaes-ext/tree/dockerfiles-13 ;\
    # git clone https://github.com/IDAES/idaes-ext.git; \
    cd idaes-ext

##
# Build all the things
##
RUN pwd; cd repo/idaes-ext; \
    bash scripts/compile_solvers.sh; \
    bash scripts/compile_libs.sh centos7