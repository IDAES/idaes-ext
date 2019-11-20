#/bin/sh

# This could be improved a lot.  We should at least split up the solvers and the AMPL user functions.
# At this point keeping it together is a bit easier to get the initial builds going.

# Make a directory to work in
export IDAES_SRC=`pwd`

# Set a few basic things
 
export CLONE_IDAES="git clone https://github.com/idaes/idaes-dev"

# get stuff
wget https://www.coin-or.org/download/source/Ipopt/$IPOPT_VER.tgz
eval $CLONE_IDAES
tar -zxvf solvers.tgz

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

# Collect files

cd $IDAES_SRC
mkdir dist
cd dist
cp ../idaes-dev/src/dist/*.so ./
cp ../license.txt ./

if [ "$(expr substr $(uname -s) 1 7)" == "MINGW64" ]; then
    # Winodws MinGW linked libraries 
    cp /mingw64/bin/libstdc++-6.dll ./
    cp /mingw64/bin/libgcc_s_seh-1.dll ./
    cp /mingw64/bin/libwinpthread-1.dll ./
    cp /mingw64/bin/libgfortran-5.dll ./
    cp /mingw64/bin/libquadmath-0.dll ./
fi

# here you pack files
# You still have the manually rename these with 
# *-{windows, linux, or darwin}-{32 or 64}.tar.gz
# I'm not building 32 bit, but who knows ma need it.
if [ "$(expr substr $(uname -s) 1 5)" == "MINGW" ]; then
  tar -czvf idaes-lib.tar.gz *.so *.dll *.txt
else
  tar -czvf idaes-lib.tar.gz *.so *.txt
fi
