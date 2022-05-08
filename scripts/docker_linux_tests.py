flavor=$1
image="idaes-ext-${flavor}-test:latest"

docker run --name test -d -i ${image} /bin/bash

docker exec test /bin/bash -c 'cd repo
    eval "$(/root/miniconda/bin/conda shell.bash hook)"
    git clone https://github.com/eslickj/idaes-pse.git
    git checkout newbin
    cd idaes-pse
    conda create -n idaes python=3.9 pip psutil
    conda activate idaes
    pip install -e .
    idaes get-extensions --extra petsc --url https://github.com/IDAES/idaes-ext/releases/download/test-release/'

docker exec test /bin/bash -c '
    cd /repo/idaes-pse/idaes
    eval "$(/root/miniconda/bin/conda shell.bash hook)"
    conda activate idaes
    pytest -m "not integration" --ignore=dmf --ignore=commands'

docker stop test
