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
cp ../license.txt ./license_lib.txt
cp ../version.txt ./version_lib.txt
sed s/"(DATE)"/`date +%Y%m%d-%H%M`/g -i version_lib.txt

# here you pack files
tar -czvf idaes-lib-${osname}-64.tar.gz *
