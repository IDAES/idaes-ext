FROM centos:centos8
MAINTAINER The CentOS Project <cloud-ops@centos.org>

RUN yum -y update; yum clean all
RUN yum -y install git make gcc gcc-c++ gcc-gfortran boost-devel cmake
RUN yum -y install wget unzip patch which file
RUN yum -y install dnf-plugins-core
RUN yum config-manager --set-enabled PowerTools
RUN dnf --enablerepo=PowerTools -y install lapack-devel
RUN mkdir repo

ENTRYPOINT ["/bin/bash"]
