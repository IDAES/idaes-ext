FROM centos:centos7
MAINTAINER The CentOS Project <cloud-ops@centos.org>

RUN yum -y update; yum clean all
RUN yum -y install git make gcc gcc-c++ gcc-gfortran boost-devel
RUN yum -y install wget unzip patch which file lapack-devel blas-devel
RUN yum -y install epel-release
RUN yum -y install cmake3
RUN ln -s /usr/bin/cmake3 /usr/bin/cmake

##
# Pull down the repo
##
RUN mkdir repo
