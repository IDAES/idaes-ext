#!/bin/sh
# First argument is OS name provided by user when running this
osname=$1; shift
if [ -z $osname ]
then
  echo "Must spcify plaform in {windows, darwin, centos6, centos7, centos8, "
  echo "  ubuntu1804, ubuntu1910, ubuntu2004}."
  exit 1
fi

# Get path of directory we're working in
export IDAES_EXT=`pwd`

# Set a few basic things
export IPOPT_BRANCH="idaes-3.13"
export IPOPT_REPO="https://github.com/idaes/Ipopt"
export PYNU_BRANCH="master"
export PYNU_REPO="https://github.com/pyomo/pyomo"
export K_AUG_BRANCH="default"
export K_AUG_REPO="https://github.com/dthierry/k_aug"
export GCC="gcc"

# Work-around for mumps gcc v10 gfortran bug
export GCC_VERSION=`$GCC -dumpversion`
if [ "$(expr substr "$GCC_VERSION" 1 2)" = "10" ]; then
  export FCFLAGS="-w -fallow-argument-mismatch -O2"
  export FFLAGS="-w -fallow-argument-mismatch -O2"
fi

mkdir coinbrew
cd coinbrew
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew

# Fetch coin-or stuff and dependencies
bash coinbrew fetch Clp --no-prompt
bash coinbrew fetch Cbc --no-prompt
bash coinbrew fetch Bonmin --no-prompt
bash coinbrew fetch Couenne --no-prompt
rm -rf Ipopt # Remove the version of Ipopt gotten as a dependency
#bash coinbrew fetch MibS@stable/1.1 --no-prompt
rm -rf ThirdParty/ASL # Remove ASL and let Ipopt have what it wants
bash coinbrew fetch $IPOPT_REPO@$IPOPT_BRANCH --no-prompt

# I let fetch get dependencies, but delete ones where I don't have permission
# to distribute, have incompatible licences, or I want a differnt version. This
# isn't really needed, but it helps keep track of what I want and what I don't
rm -rf ThirdParty/Blas
rm -rf ThirdParty/Lapack
rm -rf ThirdParty/FilterSQP
rm -rf ThirdParty/SCIP
rm -rf ThirdParty/SoPlex
rm -rf ThirdParty/glpk

# If we have the HSL stuff copy and extract it in the right place
if [ -f $IDAES_EXT/../coinhsl.zip ]
then
  # if the HSL source zip is in place...
  mkdir ThirdParty/HSL/coinhsl
  cp $IDAES_EXT/../coinhsl.zip ThirdParty/HSL/coinhsl/
  cd ThirdParty/HSL/coinhsl
  unzip coinhsl.zip
  cd $IDAES_EXT/coinbrew
  with_hsl="YES"
else
  # If the HSL isn't there, build without it.
  echo "HSL Not Available, BUILDING SOLVERS WITHOUT HSL" >&2
  with_hsl="NO"
fi

# Set PKG_CONFIG_PATH so configure scripts can find the stuff we build
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:${IDAES_EXT}/coinbrew/dist/lib/pkgconfig

echo "#########################################################################"
echo "# Thirdparty/ASL                                                        #"
echo "#########################################################################"
cd ThirdParty/ASL
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Thirdparty/Metis                                                      #"
echo "#########################################################################"
cd ThirdParty/Metis
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Thirdparty/HSL                                                        #"
echo "#########################################################################"
cd ThirdParty/HSL
./configure --disable-shared --enable-static --with-metis --prefix=$IDAES_EXT/coinbrew/dist FFLAGS="-fPIC" CFLAGS="-fPIC" CXXFLAGS="-fPIC"
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Thirdparty/Mumps                                                      #"
echo "#########################################################################"
cd ThirdParty/Mumps
./configure --disable-shared --enable-static --with-metis --prefix=$IDAES_EXT/coinbrew/dist FFLAGS="-fPIC" CFLAGS="-fPIC" CXXFLAGS="-fPIC"
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Data/Netlib                                                           #"
echo "#########################################################################"
cd Data/Netlib
./configure --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Data/Sample                                                           #"
echo "#########################################################################"
cd Data/Sample
./configure --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Data/miplib3                                                          #"
echo "#########################################################################"
cd Data/miplib3
./configure --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Ipopt ampl executables                                                #"
echo "#########################################################################"
cd Ipopt
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# CoinUtils                                                             #"
echo "#########################################################################"
cd CoinUtils
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Osi                                                                   #"
echo "#########################################################################"
cd Osi
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Clp                                                                   #"
echo "#########################################################################"
cd Clp
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Cgl                                                                   #"
echo "#########################################################################"
cd Cgl
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Cbc                                                                   #"
echo "#########################################################################"
cd Cbc
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist LDFLAGS=-fopenmp
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Bonmin                                                                #"
echo "#########################################################################"
cd Bonmin
# Two replacemnts below are a temporary fix until I can nail the problem down better
# There is a problem linking the exception class TMINLP_INVALID so replaced
# it with a similar one from IPOPT that works
sed s/"TMINLP_INVALID"/"INVALID_TNLP"/g Bonmin/src/Interfaces/BonTMINLP2TNLP.cpp > atmpfile
mv atmpfile Bonmin/src/Interfaces/BonTMINLP2TNLP.cpp
sed s/"TMINLP_INVALID"/"INVALID_TNLP"/g Bonmin/src/Interfaces/BonBranchingTQP.cpp > atmpfile
mv atmpfile Bonmin/src/Interfaces/BonBranchingTQP.cpp
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Couenne                                                               #"
echo "#########################################################################"
cd Couenne
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist LDFLAGS=-fopenmp
make
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Ipopt Shared Libraries                                                #"
echo "#########################################################################"
cd Ipopt
./configure --enable-shared --without-asl --prefix=$IDAES_EXT/coinbrew/dist
make
make install
cd $IDAES_EXT/coinbrew

# Copy files
cd $IDAES_EXT
mkdir dist-solvers
cd dist-solvers
# Executables
cp ../coinbrew/dist/bin/ipopt ./
cp ../coinbrew/dist/bin/ipopt_sens ./
cp ../coinbrew/dist/bin/clp ./
cp ../coinbrew/dist/bin/cbc ./
cp ../coinbrew/dist/bin/bonmin ./
cp ../coinbrew/dist/bin/couenne ./
# Windows *.DLL, be explicit so don't include anything we shouldn't
cp ../coinbrew/dist/bin/libipopt*.dll ./
cp ../coinbrew/dist/bin/libsipopt*.dll ./
# Linux *.so, be explicit so don't include anything we shouldn't
cp ../coinbrew/dist/lib/libipopt*.so ./
cp ../coinbrew/dist/lib/libsipopt*.so ./
# Run strip to remove unneeded symbols (it's okay that some files
#  aren't exe or libraries)
strip --strip-unneeded *

# Text information files include build time
cp ../license.txt ./
cp ../version.txt ./version_solvers.txt
sed s/"(DATE)"/`date +%Y%m%d-%H%M`/g version_solvers.txt > tmp
sed s/"(PLAT)"/${osname}/g tmp > tmp2
mv tmp2 version_solvers.txt

if [ ${osname} = "windows" ]
then
    # Winodws MinGW linked redistributable libraries
    cp /mingw64/bin/libstdc++-6.dll ./
    cp /mingw64/bin/libgcc_s_seh-1.dll ./
    cp /mingw64/bin/libwinpthread-1.dll ./
    cp /mingw64/bin/libgfortran-*.dll ./
    cp /mingw64/bin/libquadmath-0.dll ./
    cp /mingw64/bin/libgomp-1.dll ./
    cp /mingw64/bin/liblapack.dll ./
    cp /mingw64/bin/libblas.dll ./
    cp /mingw64/bin/libbz2-*.dll ./
    cp /mingw64/bin/zlib*.dll ./
fi

# Compile Pynumero

cd $IDAES_EXT
git clone $PYNU_REPO
cd pyomo
git checkout $PYNU_BRANCH
cd pyomo/contrib/pynumero/src
mkdir build
cd build
if [ "$(expr substr $(uname -s) 1 7)" = "MINGW64" ]
then
  cmake -DENABLE_HSL=no -DIPOPT_DIR=$IDAES_EXT/coinbrew/dist -G"MSYS Makefiles" ..
else
  cmake .. -DENABLE_HSL=no -DIPOPT_DIR=$IDAES_EXT/coinbrew/dist
fi
make
cp libpynumero_ASL* $IDAES_EXT/dist-solvers

# Compile k_aug
cd $IDAES_EXT
git clone $K_AUG_REPO
cp ./scripts/k_aug_CMakeLists.txt ./k_aug/CMakeLists.txt
cd k_aug
git checkout $K_AUG_BRANCH
if [ "$(expr substr $(uname -s) 1 7)" = "MINGW64" ]
then
  cmake -DWITH_MINGW=ON -DCMAKE_C_COMPILER=gcc -G"MSYS Makefiles" .
else
  cmake -DCMAKE_C_COMPILER=gcc .
fi
make
cp bin/k_aug* $IDAES_EXT/dist-solvers
cp dot_sens* $IDAES_EXT/dist-solvers

# here you pack files
cd $IDAES_EXT/dist-solvers
tar -czvf idaes-solvers-${osname}-64.tar.gz *
echo "HSL Present: ${with_hsl}"
