FROM eslickj/idaes-ext-ubuntu1804-build:latest

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
    bash scripts/compile_solvers.sh ubuntu1804; \
    bash scripts/compile_libs.sh ubuntu1804
