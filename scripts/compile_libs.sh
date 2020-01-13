#/bin/sh

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
mkdir dist-lib
cd dist-lib
cp ../src/dist/*.so ./
cp ../license.txt ./

if [ "$(expr substr $(uname -s) 1 7)" == "MINGW64" ]; then
    # Winodws MinGW linked libraries
    cp /mingw64/bin/libstdc++-6.dll ./
    cp /mingw64/bin/libgcc_s_seh-1.dll ./
    cp /mingw64/bin/libwinpthread-1.dll ./
    cp /mingw64/bin/libgfortran-5.dll ./
    cp /mingw64/bin/libquadmath-0.dll ./
    cp /mingw64/bin/liblapack.dll ./
    cp /mingw64/bin/libblas.dll ./
fi

# here you pack files
tar -czvf idaes-lib.tar.gz *
