#/bin/sh
# First argument is OS name
osname=$1; shift

# Make a directory to work in
export IDAES_EXT=`pwd`

# get stuff
wget https://ampl.com/netlib/ampl/solvers.tgz
tar -zxvf solvers.tgz

# Compile ASL, warnings about files existing seem to be okay

cd solvers
./configure
make
export ASL_BUILD=`pwd`/sys.`uname -m`.`uname -s`
cd $IDAES_EXT

# Compile IDAES function libraries

cd src
make
cd $IDAES_EXT

# Collect files

cd $IDAES_EXT
rm -rf ./dist-lib
mkdir dist-lib
cd dist-lib
cp ../src/dist/*.so ./
cp ../license.txt ./

if [ "$(expr substr $(uname -s) 1 7)" == "MINGW64" ]
then
    # Winodws MinGW linked libraries
    cp /mingw64/bin/libstdc++-6.dll ./
    cp /mingw64/bin/libgcc_s_seh-1.dll ./
    cp /mingw64/bin/libwinpthread-1.dll ./
    cp /mingw64/bin/libgfortran-4.dll ./
    cp /mingw64/bin/libquadmath-0.dll ./
    cp /mingw64/bin/libgomp-1.dll ./
    cp /mingw64/bin/liblapack.dll ./
    cp /mingw64/bin/libblas.dll ./
fi

# Compile Pynumero

cd $IDAES_EXT
git clone https://github.com/idaes/pyomo pyomo
cd pyomo
git checkout IDAES
cd pyomo/contrib/pynumero/cmake/third_party/ASL
sh ./getASL.sh
cd solvers
sh ./configurehere
sed -e s/-DNo_dtoa//g -i Makefile
make
cd ../../../
mkdir build
cd build
if [ "$(expr substr $(uname -s) 1 7)" == "MINGW64" ]
then
  cmake -G"MSYS Makefiles" ..
else
  cmake ..
fi
make
cp asl_interface/libpynumero_ASL* $IDAES_EXT/dist-lib
#cp sparse_utils/libpynumero_SPARSE* $IDAES_EXT/dist-lib
cd $IDAES_EXT/dist-lib

# here you pack files
tar -czvf idaes-lib-${osname}-64.tar.gz *
