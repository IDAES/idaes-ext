FROM centos:centos7
MAINTAINER The CentOS Project <cloud-ops@centos.org>

ENV repo https://github.com/eslickj/idaes-ext.git
ENV branch dock

RUN yum -y update; yum clean all
RUN yum -y install git make gcc gcc-c++ gcc-gfortran boost-devel
RUN yum -y install wget

##
# Pull down the repo
##
RUN mkdir repo

COPY ./extras/ /repo/

RUN cd repo; \
    git clone ${repo} ; \
    cd idaes-ext ; \
    git checkout ${branch}; \
    echo "cloned repo"

##
# Build all the things
##
RUN pwd; cd repo/idaes-ext; \
    # don't need to compile solvers for now
    bash scripts/compile_libs.sh centos7; \
    bash scripts/compile_solvers.sh centos7
