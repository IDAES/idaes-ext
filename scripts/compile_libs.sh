#!/bin/sh
# First argument is OS name
osname=$1; shift
if [ -z $osname ]
then
  echo "Must specify platform in {windows, darwin, el8, el9, ubuntu2004, ubuntu2204, ubuntu2404}."
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

if [ ${osname} = "darwin" ]; then
  # Homebrew boost is required
  export BOOST=$($HOMEBREW_PREFIX/bin/brew --prefix boost)/include
fi

make BOOST=$BOOST
cd $IDAES_EXT

# Collect files
rm -rf ./dist-functions
mkdir dist-functions dist-functions/lib
#cd dist-functions
cp ./src/dist/*.so ./dist-functions/lib/
if [ ${osname} = "windows" ]; then
  mv dist-functions/lib/functions.so dist-functions/lib/functions.dll
  mv dist-functions/lib/general_helmholtz_external.so dist-functions/lib/general_helmholtz_external.dll
  mv dist-functions/lib/cubic_roots.so dist-functions/lib/cubic_roots.dll
fi
if [ ${osname} = "darwin" ]; then
  mv dist-functions/lib/functions.so dist-functions/lib/functions.dylib
  mv dist-functions/lib/general_helmholtz_external.so dist-functions/lib/general_helmholtz_external.dylib
  mv dist-functions/lib/cubic_roots.so dist-functions/lib/cubic_roots.dylib
fi
cp ./license.txt ./dist-functions/license_functions.txt
cp ./version.txt ./dist-functions/version_functions.txt
mkdir ./dist-functions/lib/helm_data
cp ./src/dist/param_data/*.json ./dist-functions/lib/helm_data/
cp ./src/dist/param_data/*.nl ./dist-functions/lib/helm_data/
cp ./src/dist/param_data/*.py ./dist-functions/lib/helm_data/
sed s/"(DATE)"/`date +%Y%m%d-%H%M`/g dist-functions/version_functions.txt > tmp
sed s/"(PLAT)"/${osname}-${MNAME}/g tmp > tmp2
mv tmp2 dist-functions/version_functions.txt
rm tmp

# here you pack files
cd dist-functions
tar -czvf idaes-functions-${osname}-${MNAME}.tar.gz *
cd $IDAES_EXT
