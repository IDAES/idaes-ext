#!/bin/sh
# First argument is OS name
osname=$1; shift
if [ -z $osname ]
then
  echo "Must spcify plaform in {windows, darwin, el7, el8, ubuntu1804, ubuntu2004}."
  exit 1
fi

# Get the machine type
export MNAME=`uname -m`

# Make a directory to work in
export IDAES_EXT=`pwd`

# Run this after solvers are compiled, uses the ASL header from solver build
export ASL_HEADERS=$IDAES_EXT/coinbrew/dist/include/coin-or/asl
export ASL_LIBRARIES=$IDAES_EXT/coinbrew/dist/lib

# Compile IDAES function libraries
cd $IDAES_EXT/src
if [ ${osname} = "windows" ]; then
  export PATH=$PATH:$IDAES_EXT/coinbrew/dist/lib64
fi
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$IDAES_EXT/coinbrew/dist/lib64
export IDAES_HELMHOLTZ_DATA_PATH=$IDAES_EXT/src/general_helmholtz/param_data/
export IDAES_HELMHOLTZ_TEST_DATA_PATH=$IDAES_EXT/src/general_helmholtz/test_data/

make
cd $IDAES_EXT

# Collect files
rm -rf ./dist-lib
mkdir dist-lib
cd dist-lib
cp ../src/dist/*.so ./
if [ ${osname} = "windows" ]; then
  mv functions.so functions.dll
  mv general_helmholtz_external.so general_helmholtz_external.dll
  mv cubic_roots.so cubic_roots.dll
fi
if [ ${osname} = "darwin" ]; then
  mv functions.so functions.dylib
  mv general_helmholtz_external.so general_helmholtz_external.dylib
  mv cubic_roots.so cubic_roots.dylib
fi
cp ../license.txt ./license_lib.txt
cp ../version.txt ./version_lib.txt
mkdir ./helm_data
cp ../src/dist/param_data/*.json ./helm_data/
cp ../src/dist/param_data/*.nl ./helm_data/
cp ../src/dist/param_data/*.py ./helm_data/
sed s/"(DATE)"/`date +%Y%m%d-%H%M`/g version_lib.txt > tmp
sed s/"(PLAT)"/${osname}-${MNAME}/g tmp > tmp2
mv tmp2 version_lib.txt
rm tmp

# here you pack files
tar -czvf idaes-lib-${osname}-${MNAME}.tar.gz *
