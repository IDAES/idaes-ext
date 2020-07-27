FROM eslickj/idaes-ext-ubuntu1910-test:latest

ARG repo=https://github.com/idaes/idaes-pse.git
ARG branch=master

RUN cd repo; \
    eval "$(/root/miniconda/bin/conda shell.bash hook)"; \
    git clone ${repo}; \
    cd idaes-pse; \
    git checkout ${branch}; \
    echo "cloned repo"; \
    conda create -n idaes python=3.7 pip psutil; \
    conda activate idaes; \
    pip install -e .; \
    cd ..; \
    idaes get-extensions --url https://github.com/IDAES/idaes-ext/releases/download/test-release/ --platform ubuntu1910

RUN cd repo; \
    eval "$(/root/miniconda/bin/conda shell.bash hook)"; \
    conda activate idaes; \
    cd idaes-pse/idaes; \
    pytest tests/test_have_pynumero.py; \
    pytest tests/test_ipopt.py;\
    cd ./generic_models; \
    pytest; \
    cd ../core; \
    pytest; \
    cd ../power_generation; \
    pytest
