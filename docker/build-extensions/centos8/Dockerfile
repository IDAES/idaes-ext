FROM eslickj/idaes-ext-centos8-build:latest

ARG repo=https://github.com/idaes/idaes-ext.git
ARG branch=master

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
    bash scripts/compile_solvers.sh centos8; \
    bash scripts/compile_libs.sh centos8
