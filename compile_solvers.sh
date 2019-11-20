#!/bin/sh

# This could be improved a lot.  We should at least split up the solvers and the AMPL user functions.
# At this point keeping it together is a bit easier to get the initial builds going.

# Make a directory to work in
export IDAES_SRC=`pwd`

# Set a few basic things
 
export CLONE_IPOPT="git clone --single-branch --branch idaes-3.12.13 https://github.com/idaes/Ipopt"
export IPOPT_VER=Ipopt-3.12.13
export PATCH_IPOPT="cp Ipopt/Ipopt/src/Algorithm/IpIpoptAlg.cpp Ipopt-3.12.13/Ipopt/src/Algorithm/IpIpoptAlg.cpp"

# get stuff
wget https://www.coin-or.org/download/source/Ipopt/$IPOPT_VER.tgz
eval $CLONE_IPOPT
tar -zxvf $IPOPT_VER.tgz
eval $PATCH_IPOPT

# Compile IPOPT

cd $IPOPT_VER/ThirdParty
cd ./Lapack
./get.Lapack
cd ../Blas
./get.Blas
cd ../ASL
./get.ASL
cd ../Mumps
./get.Mumps
cd ../Metis
./get.Metis
cd ../HSL
cp -r $IDAES_SRC/coinhsl ./
cd ../../
./configure --enable-shared=no
make

# Collect files

cd $IDAES_SRC
mkdir dist-solvers
cd dist-solvers
cp ../Ipopt-3.12.13/Ipopt/src/Apps/AmplSolver/ipopt ./
cp ../license.txt ./

# here you pack files
tar -czvf idaes-solvers.tar.gz *
