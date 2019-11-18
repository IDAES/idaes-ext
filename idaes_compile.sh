#!/bin/sh

# This could be improved a lot.  We should at least split up the solvers and the AMPL user functions.
# At this point keeping it together is a bit easier to get the initial builds going.

# Make a directory to work in
export IDAES_SRC=`pwd`

# Set a few basic things
 
export CLONE_IDAES="git clone https://github.com/idaes/idaes-dev"
export CLONE_IPOPT="git clone --single-branch --branch idaes-3.12.13 https://github.com/idaes/Ipopt"
export IPOPT_VER=Ipopt-3.12.13
export PATCH_IPOPT="cp Ipopt/Ipopt/src/Algorithm/IpIpoptAlg.cpp Ipopt-3.12.13/Ipopt/src/Algorithm/IpIpoptAlg.cpp"

# get stuff
wget https://ampl.com/netlib/ampl/solvers.tgz
wget https://www.coin-or.org/download/source/Ipopt/$IPOPT_VER.tgz
eval $CLONE_IDAES
eval $CLONE_IPOPT
tar -zxvf solvers.tgz
tar -zxvf $IPOPT_VER.tgz
eval $PATCH_IPOPT

# Compile ASL, warnings about files existing seem to be okay

cd solvers
./configure
make
export ASL_BUILD=`pwd`/sys.`uname -m`.`uname -s`
cd $IDAES_SRC

# Compile IDAES function libraries

cd idaes-dev/src
make
cd $IDAES_SRC

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
mkdir dist
cd dist
cp ../idaes-dev/src/dist/*.so ./
cp ../Ipopt-3.12.13/Ipopt/src/Apps/AmplSolver/ipopt ./
cp $IDAES_SRC/license.txt ./

# Winodws MinGW linked libraries 
# cp /mingw64/bin/libstdc++-6.dll ./
# cp /mingw64/bin/libgcc_s_seh-1.dll ./
# cp /mingw64/bin/libwinpthread-1.dll ./
# cp /mingw64/bin/libgfortran-5.dll ./
# cp /mingw64/bin/libquadmath-0.dll ./

# here you zip files

zip idaes-lib.zip *.so *.dll
zip idaes-solvers.zip ipopt*
