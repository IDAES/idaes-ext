FROM eslickj/idaes-ext-windows-test:latest

ARG repo=https://github.com/idaes/idaes-pse.git
ARG branch=master

SHELL ["powershell", "-command"]

WORKDIR c:/windows/temp

RUN cd \repo; \
    conda activate idaes; \
    git clone $Env:repo;\
    cd idaes-pse; \
    git checkout $env:branch; \
    pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org -r requirements.txt; \
    pip install -e .; \
    idaes get-extensions --url https://github.com/IDAES/idaes-ext/releases/download/test-release/

RUN cd \repo; \
    conda activate idaes; \
    cd idaes; \
    pytest tests\test_have_pynumero.py; \
    pytest tests\test_ipopt.py;\
    cd generic_models; \
    pytest; \
    cd ../core; \
    pytest
    #cd ../power_generation; \
    #pytest
