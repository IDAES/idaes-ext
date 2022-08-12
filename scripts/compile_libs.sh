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
export ASL_BUILD=$IDAES_EXT/coinbrew/dist/include/coin-or/asl
export ADOLC_INC=$IDAES_EXT/coinbrew/dist/include
export ADOLC_LIB=$IDAES_EXT/coinbrew/dist/lib64

# Compile IDAES function libraries
cd $IDAES_EXT/src
if [ ${osname} = "darwin" ]; then
  export BOOST_HEADER=/opt/homebrew/include
fi
if [ ${osname} = "darwin" ]; then
  export PATH=$PATH:$IDAES_EXT/coinbrew/dist/lib64
fi
export LD_LIBRARY_PATH=$IDAES_EXT/coinbrew/dist/lib64
make
cd $IDAES_EXT

# Collect files
rm -rf ./dist-lib
mkdir dist-lib
cd dist-lib
cp ../src/dist/*.so ./
cp ../license.txt ./license_lib.txt
cp ../version.txt ./version_lib.txt
sed s/"(DATE)"/`date +%Y%m%d-%H%M`/g version_lib.txt > tmp
sed s/"(PLAT)"/${osname}-${MNAME}/g tmp > tmp2
mv tmp2 version_lib.txt
rm tmp

# here you pack files
tar -czvf idaes-lib-${osname}-${MNAME}.tar.gz *
