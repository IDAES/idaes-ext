#!/bin/sh

# Make a directory to work in
export IDAES_SRC=`pwd`

# Set a few basic things

export IPOPT_BRANCH="idaes-3.12.13"
export IPOPT_REPO="https://github.com/idaes/Ipopt"

mkdir coinbrew
cd coinbrew
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
bash coinbrew fetch $IPOPT_REPO:$IPOPT_BRANCH --no-prompt
echo -n >ThirdParty/HSL/.build
cp -r $IDAES_SRC/coinhsl ThirdParty/HSL/coinhsl
bash coinbrew build $IPOPT_REPO:$IPOPT_BRANCH --no-prompt

cd $IDAES_SRC
mkdir dist-solvers
cd dist-solvers
cp ../coinbrew/dist/bin/ipopt ./
cp ../license.txt ./

# here you pack files
tar -czvf idaes-solvers.tar.gz *
