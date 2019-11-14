
# Make a directory to work in

cd $IDAES_SRC

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
# Waiting to resolve METIS license issue
#cd ../Metis
#./get.Metis
cd ../HSL
cp -r $IDAES_SRC/coinhsl ./
cd ../../
./configure CFLAGS="-static"
make LDFLAGS="-all_static"

# Collect files

cd $IDAES_SRC
mkdir dist
cd dist
cp ../idaes-dev/src/dist/*.so ./
cp ../Ipopt-3.12.13/Ipopt/src/Apps/AmplSolver/ipopt ./

# here you zip files

zip idaes-bin.zip *
