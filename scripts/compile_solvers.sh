#!/bin/sh

# Set a working directory (I'm assuming this will be the idaes-ext repo dir)
export IDAES_EXT=`pwd`

# Set a few basic things

export IPOPT_BRANCH="idaes-3.13.0"
export IPOPT_REPO="https://github.com/idaes/Ipopt"


mkdir coinbrew
cd coinbrew
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
bash coinbrew fetch $IPOPT_REPO:$IPOPT_BRANCH --no-prompt
echo -n >ThirdParty/HSL/.build
cp -r $IDAES_EXT/coinhsl ThirdParty/HSL/coinhsl
bash coinbrew build $IPOPT_REPO:$IPOPT_BRANCH --no-prompt

# I think there is an easier way, but I'm going to explicitly do these one-by-one
bash coinbrew fetch Clp:stable/1.17 --no-prompt
bash coinbrew build Clp --no-prompt

bash coinbrew fetch Cbc:stable/2.10 --no-prompt
bash coinbrew build Cbc --no-prompt

bash coinbrew fetch Bonmin:stable/1.8 --no-prompt
bash coinbrew build Bonmin --no-prompt

bash coinbrew fetch Couenne:stable/0.5 --no-prompt
bash coinbrew build Couenne --no-prompt


cd $IDAES_EXT
mkdir dist-solvers
cd dist-solvers
cp ../coinbrew/dist/bin/ipopt ./
cp ../coinbrew/dist/bin/clp ./
cp ../coinbrew/dist/bin/ipopt_sens ./
cp ../license.txt ./

# here you pack files
tar -czvf idaes-solvers.tar.gz *
