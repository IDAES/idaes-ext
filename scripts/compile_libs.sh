#/bin/sh
# First argument is OS name
osname=$1; shift

# Make a directory to work in
export IDAES_EXT=`pwd`

# get stuff
wget https://coin-or-tools.github.io/ThirdParty-ASL/solvers-20180528.tgz
mv solvers-20180528.tgz solvers.tgz
tar -zxvf solvers.tgz

# Compile ASL, warnings about files existing seem to be okay

cd solvers

# adowling laptop
# ./configure CC=gcc-10 F77=gfortran-10

# adowling desktop
./configure CC=gcc-9 F77=gfortran-9


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
cp ../license.txt ./license_lib.txt
cp ../version.txt ./version_lib.txt
sed s/"(DATE)"/`date +%Y%m%d-%H%M`/g version_lib.txt > tmp
sed s/"(PLAT)"/${osname}/g tmp > tmp2
mv tmp2 version_lib.txt

# here you pack files
tar -czvf idaes-lib-${osname}-64.tar.gz *
