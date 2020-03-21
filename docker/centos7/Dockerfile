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
    # debug mode: clone from a fork
    # git clone https://github.com/dangunter/idaes-ext.git ; \
    git clone https://github.com/eslickj/idaes-ext.git ; \
    # normal mode: clone from IDAES
    # git clone https://github.com/IDAES/idaes-ext.git; \
    cd idaes-ext ; \
    # debug mode: check out PR branch
    # git checkout dockerfiles-13; \
    echo "cloned repo"

##
# Build all the things
##
RUN pwd; cd repo/idaes-ext; \
    bash scripts/compile_solvers.sh centos7; \
    bash scripts/compile_libs.sh centos7
