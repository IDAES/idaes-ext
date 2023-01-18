#!/bin/sh
# First argument is OS name provided by user when running this
osname=$1;
if [ -z $osname ]
then
  echo "Must spcify plaform in {windows, darwin, el7, el8, "
  echo "  ubuntu180, ubuntu2004, ubuntu2204}."
  exit 1
fi

# Get the machine type
export MNAME=`uname -m`

# Get path of directory we're working in
export IDAES_EXT=`pwd`

if [ -f $IDAES_EXT/../coinhsl.zip ]
then
  echo "HSL: YES"
else
  echo "HSL: NO"
fi

# Set a few basic things
export IPOPT_BRANCH="idaes-3.13"
export IPOPT_REPO="https://github.com/idaes/Ipopt"
export IPOPT_L1_BRANCH="restoration_mod"
export IPOPT_L1_REPO="https://github.com/idaes/Ipopt"
export PYNU_BRANCH="main"
export PYNU_REPO="https://github.com/pyomo/pyomo"
export K_AUG_BRANCH="default"
export K_AUG_REPO="https://github.com/dthierry/k_aug"
export CC="gcc"
export CXX="g++"

# set PETSc location
if [ ${osname} = "windows" ]
then
  export PETSC_DIR=/c/repo/petsc-dist
elif [ ${osname} = "darwin" ]; then
  export PETSC_DIR="$HOME/src/petsc-dist"
else
  export PETSC_DIR=/repo/petsc-dist
fi
export PETSC_ARCH=""

mkdir coinbrew
cd coinbrew

if [ ${osname} = "darwin" ]; then
  curl --output coinbrew https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
  export CC="gcc-12"
  export CXX="g++-12"
else
  wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
fi

# Work-around for mumps gcc v10 gfortran bug
GFORT_VERSION=`gfortran -dumpversion`
GFMV=(${GFORT_VERSION//./ })
if [ ${GFMV[0]} -ge 10 ]; then
  export FCFLAGS="-w -fallow-argument-mismatch -O2"
  export FFLAGS="-w -fallow-argument-mismatch -O2"
fi

# Fetch coin-or stuff and dependencies
SKIP_PKGS='ThirdParty/Lapack ThirdParty/Blas ThirdParty/glpk ThirdParty/Metis ThirdParty/Mumps'
bash coinbrew fetch Clp --no-prompt --skip "$SKIP_PKGS"
bash coinbrew fetch Cbc --no-prompt --skip "$SKIP_PKGS"
SKIP_PKGS="$SKIP_PKGS Cbc Clp Cgl Osi"
bash coinbrew fetch Bonmin --no-prompt --skip "$SKIP_PKGS"
bash coinbrew fetch Couenne --no-prompt --skip "$SKIP_PKGS"
# Patch Couenne to fix: error: static assertion failed: comparison object must be invocable as const
cd Couenne
cp $IDAES_EXT/scripts/CouenneMatrix.hpp.patch ./
cp $IDAES_EXT/scripts/CouenneProblem.hpp.patch ./
patch Couenne/src/problem/CouenneProblem.hpp < CouenneProblem.hpp.patch
patch Couenne/src/cut/sdpcuts/CouenneMatrix.hpp < CouenneMatrix.hpp.patch
cd ..
rm -rf Ipopt # Remove the version of Ipopt gotten as a dependency
bash coinbrew fetch $IPOPT_L1_REPO@$IPOPT_L1_BRANCH --no-prompt --skip "$SKIP_PKGS"
mv ./Ipopt ./Ipopt_l1
rm -rf ThirdParty/ASL # Remove ASL and let Ipopt have what it wants
if [ ${osname} = "el7" ]; then 
  # Seems now the git autostash option is used, but not in the older git in el7. To prevent failure to get Mumps, just delete Thirdparty
  # Looks like the only things in therd party are thing that Ipopt gets, so this should be ok.
  rm -rf ./Thirdparty/*
fi
bash coinbrew fetch $IPOPT_REPO@$IPOPT_BRANCH --no-prompt --skip 'ThirdParty/Lapack ThirdParty/Blas ThirdParty/glpk'
cp -r Ipopt Ipopt_share

# Make sure I don't include any dependencies I don't want
rm -rf ThirdParty/Blas
rm -rf ThirdParty/Lapack
rm -rf ThirdParty/FilterSQP
rm -rf ThirdParty/SCIP
rm -rf ThirdParty/SoPlex
rm -rf ThirdParty/Glpk

if [ -f $IDAES_EXT/../metis-4.0-novariadic.tar.gz ]; then
  rm -rf ./ThirdParty/Metis/metis-4.0/*
  cp $IDAES_EXT/../metis-4.0-novariadic.tar.gz ./ThirdParty/Metis/metis-4.0/
  cd ThirdParty/Metis/metis-4.0
  tar -zxvf metis-4.0-novariadic.tar.gz
  rm metis-4.0-novariadic.tar.gz
  echo "#########################################################################"
  echo "# Use novariadic metis 4.0.3                                            #"
  echo "#########################################################################"
  make $PARALLEL
  cd $IDAES_EXT/coinbrew
fi

######
# Set PKG_CONFIG_PATH so configure scripts can find the stuff we build
#######
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:${IDAES_EXT}/coinbrew/dist/lib/pkgconfig

echo "#########################################################################"
echo "# Get coinhsl.zip if available                                          #"
echo "#########################################################################"
# If we have the HSL stuff copy and extract it in the right place
if [ -f $IDAES_EXT/../coinhsl.zip ]; then
  # if the HSL source zip is in place...
  mkdir ThirdParty/HSL/coinhsl
  cp $IDAES_EXT/../coinhsl.zip ThirdParty/HSL/coinhsl/
  cd ThirdParty/HSL/coinhsl
  unzip coinhsl.zip
  cd $IDAES_EXT/coinbrew
  echo "HSL is available, building with HSL"
  with_hsl="YES"
else
  # If the HSL isn't there, build without it.
  echo "HSL Not Available, BUILDING SOLVERS WITHOUT HSL" >&2
  with_hsl="NO"
fi

echo "#########################################################################"
echo "# Thirdparty/ASL                                                        #"
echo "#########################################################################"
cd ThirdParty/ASL
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Thirdparty/Metis                                                      #"
echo "#########################################################################"
cd ThirdParty/Metis
./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist \
  --prefix=$IDAES_EXT/coinbrew/dist FFLAGS="-fPIC" CFLAGS="-fPIC" CXXFLAGS="-fPIC"
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Thirdparty/HSL                                                        #"
echo "#########################################################################"
cd ThirdParty/HSL
./configure --disable-shared --enable-static --with-metis \
  --prefix=$IDAES_EXT/coinbrew/dist FFLAGS="-fPIC" CFLAGS="-fPIC" CXXFLAGS="-fPIC"
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Thirdparty/Mumps                                                      #"
echo "#########################################################################"
cd ThirdParty/Mumps
./configure --disable-shared --enable-static --with-metis \
 --prefix=$IDAES_EXT/coinbrew/dist FFLAGS="-fPIC" CFLAGS="-fPIC" CXXFLAGS="-fPIC"
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Ipopt ampl executables                                                #"
echo "#########################################################################"
cd Ipopt
./configure --disable-shared --enable-static --with-mumps --with-hsl \
  --prefix=$IDAES_EXT/coinbrew/dist
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Ipopt_L1 ampl executables                                             #"
echo "#########################################################################"
cd Ipopt_l1
if [ ${osname} = "el7" ]; then
  ./configure --disable-shared --enable-static --with-mumps --with-hsl \
    ADD_CXXFLAGS="-std=c++11" \
    --prefix=$IDAES_EXT/coinbrew/dist_l1
else
  ./configure --disable-shared --enable-static --with-mumps --with-hsl \
    --prefix=$IDAES_EXT/coinbrew/dist_l1
fi
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# CoinUtils                                                             #"
echo "#########################################################################"
cd CoinUtils
if [ "$MNAME" = "aarch64" ]; then
  # MNAME of darwin is arm64, so this is linux only
  ./configure --build=aarch64-unknown-linux-gnu --disable-shared --enable-static \
    --prefix=$IDAES_EXT/coinbrew/dist
else
  ./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
fi
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Osi                                                                   #"
echo "#########################################################################"
cd Osi
if [ "$MNAME" = "aarch64" ]; then
  # MNAME of darwin is arm64, so this is linux only
  ./configure --build=aarch64-unknown-linux-gnu --disable-shared --enable-static \
    --prefix=$IDAES_EXT/coinbrew/dist
else
  ./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
fi
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Clp                                                                   #"
echo "#########################################################################"
cd Clp
if [ "$MNAME" = "aarch64" ]; then
  # MNAME of darwin is arm64, so this is linux only
  ./configure --build=aarch64-unknown-linux-gnu --disable-shared --enable-static \
   --prefix=$IDAES_EXT/coinbrew/dist
else
  ./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
fi
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Cgl                                                                   #"
echo "#########################################################################"
cd Cgl
if [ "$MNAME" = "aarch64" ]; then
  # MNAME of darwin is arm64, so this is linux only
  ./configure --build=aarch64-unknown-linux-gnu --disable-shared --enable-static \
    --prefix=$IDAES_EXT/coinbrew/dist
else
  ./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
fi
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Cbc                                                                   #"
echo "#########################################################################"
cd Cbc
if [ "$MNAME" = "aarch64" ]; then
  # MNAME of darwin is arm64, so this is linux only
  ./configure --build=aarch64-unknown-linux-gnu --disable-shared --enable-static \
    --prefix=$IDAES_EXT/coinbrew/dist
else
  ./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist
fi
make $PARALLEL
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
if [ "$MNAME" = "aarch64" ]; then
  # MNAME of darwin is arm64, so this is linux only
  ./configure --build=aarch64-unknown-linux-gnu --disable-shared --enable-static \
    --prefix=$IDAES_EXT/coinbrew/dist LDFLAGS=-fopenmp
else
  ./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist \
    LDFLAGS=-fopenmp
fi
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Couenne                                                               #"
echo "#########################################################################"
cd Couenne
if [ "$MNAME" = "aarch64" ]; then
  # MNAME of darwin is arm64, so this is linux only
  ./configure --build=aarch64-unknown-linux-gnu --disable-shared --enable-static \
    --prefix=$IDAES_EXT/coinbrew/dist LDFLAGS=-fopenmp
else
  ./configure --disable-shared --enable-static --prefix=$IDAES_EXT/coinbrew/dist \
    LDFLAGS=-fopenmp
fi
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Ipopt Shared Libraries                                                #"
echo "#########################################################################"
cd Ipopt_share
./configure --enable-shared --disable-static --without-asl --disable-java \
  --with-mumps --with-hsl --prefix=$IDAES_EXT/coinbrew/dist-share
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Copy Coin Solver Files to dist-solvers                                #"
echo "#########################################################################"
cd $IDAES_EXT
mkdir dist-solvers
cd dist-solvers
# Executables
if [ ${osname} = "windows" ]; then
  # windows
  cp ../coinbrew/dist_l1/bin/ipopt.exe ./ipopt_l1.exe
  cp ../coinbrew/dist_l1/bin/ipopt_sens.exe ./ipopt_sens_l1.exe
  # Explicitly only get ipopt so we don't get anything we shouldn't
  cp ../coinbrew/dist-share/bin/libipopt*.dll ./
  cp ../coinbrew/dist-share/bin/libsipopt*.dll ./
elif [ ${osname} = "darwin" ]; then
  cp ../coinbrew/dist_l1/bin/ipopt ./ipopt_l1
  cp ../coinbrew/dist_l1/bin/ipopt_sens ./ipopt_sens_l1
  # Explicitly only get ipopt so we don't get anything we shouldn't
  cp ../coinbrew/dist-share/lib/libipopt*.dylib ./
  cp ../coinbrew/dist-share/lib/libsipopt*.dylib ./
else
  # linux
  cp ../coinbrew/dist_l1/bin/ipopt ./ipopt_l1
  cp ../coinbrew/dist_l1/bin/ipopt_sens ./ipopt_sens_l1
  # Explicitly only get ipopt so we don't get anything we shouldn't
  cp ../coinbrew/dist-share/lib/libipopt*.so ./
  cp ../coinbrew/dist-share/lib/libsipopt*.so ./
fi
cp ../coinbrew/dist/bin/ipopt ./
cp ../coinbrew/dist/bin/ipopt_sens ./
cp ../coinbrew/dist/bin/clp ./
cp ../coinbrew/dist/bin/cbc ./
cp ../coinbrew/dist/bin/bonmin ./
cp ../coinbrew/dist/bin/couenne ./
# Run strip to remove unneeded symbols (it's okay that some files
#  aren't exe or libraries)
strip --strip-unneeded *

cd ../
cp -r ./coinbrew/dist-share ./dist-share
cd $IDAES_EXT/dist-share
tar -czvf idaes-local-${osname}-${MNAME}.tar.gz *

cd $IDAES_EXT/dist-solvers/
cp ../dist-share/idaes-local-${osname}-${MNAME}.tar.gz ./

echo "#########################################################################"
echo "# Copy License and Version Files to dist-solvers                        #"
echo "#########################################################################"
# Text information files include build time
cp ../license.txt ./
cp ../version.txt ./version_solvers.txt
sed s/"(DATE)"/`date +%Y%m%d-%H%M`/g version_solvers.txt > tmp
sed s/"(PLAT)"/${osname}-${MNAME}/g tmp > tmp2
mv tmp2 version_solvers.txt
rm tmp

#
# Copy some linked libraries from homebrew or mingw, covered by gcc runtime
# exception https://www.gnu.org/licenses/gcc-exception-3.1.en.html license
# info is included in the license text file.
#
echo "#########################################################################"
echo "# Copy GCC/MinGW Runtime Libraries to dist-solvers                      #"
echo "#########################################################################"
if [ ${osname} = "windows" ]; then
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
    cp /mingw64/bin/libssp*.dll ./
fi

if [ ${osname} = "darwin" ]; then
    # some libraries from homebrew
    cp /opt/homebrew/opt/gcc/lib/gcc/12/libgfortran.5.dylib ./
    cp /opt/homebrew/opt/gcc/lib/gcc/12/libgcc_s.1.1.dylib ./
    cp /opt/homebrew/opt/gcc/lib/gcc/12/libstdc++.6.dylib ./
    cp /opt/homebrew/opt/gcc/lib/gcc/12/libgomp.1.dylib ./
    cp /opt/homebrew/opt/gcc/lib/gcc/12/libquadmath.0.dylib ./
fi

echo "#########################################################################"
echo "# Pynumero                                                              #"
echo "#########################################################################"

if [ ${osname} = "darwin" ]; then
  export CC="cc"
  export CXX="c++"
fi

cd $IDAES_EXT
git clone $PYNU_REPO
cd pyomo
git checkout $PYNU_BRANCH
cd pyomo/contrib/pynumero/src
mkdir build
cd build
if [ ${osname} = "windows" ]
then
  cmake -DENABLE_HSL=no -DIPOPT_DIR=$IDAES_EXT/coinbrew/dist -G"MSYS Makefiles" ..
else
  cmake .. -DENABLE_HSL=no -DIPOPT_DIR=$IDAES_EXT/coinbrew/dist
fi
make $PARALLEL
cp libpynumero_ASL* $IDAES_EXT/dist-solvers

echo "#########################################################################"
echo "# k_aug, dotsens                                                        #"
echo "#########################################################################"
cd $IDAES_EXT
git clone $K_AUG_REPO
cp ./scripts/k_aug_CMakeLists.txt ./k_aug/CMakeLists.txt
cd k_aug
git checkout $K_AUG_BRANCH
if [ ${osname} = "windows" ]
then
  cmake -DWITH_MINGW=ON -DCMAKE_C_COMPILER=$CC -G"MSYS Makefiles" .
else
  cmake -DCMAKE_C_COMPILER=$CC .
fi
make $PARALLEL
cp bin/k_aug* $IDAES_EXT/dist-solvers
cp dot_sens* $IDAES_EXT/dist-solvers

echo "#########################################################################"
echo "# PETSc                                                                 #"
echo "#########################################################################"
export ASL_INC=$IDAES_EXT/coinbrew/dist/include/coin-or/asl
export ASL_LIB=$IDAES_EXT/coinbrew/dist/lib/libcoinasl.a
cd $IDAES_EXT/petsc
make $PARALLEL
make py
mkdir $IDAES_EXT/dist-petsc
if [ ${osname} = "windows" ]
then
  cp petsc.exe $IDAES_EXT/dist-petsc
else
  cp petsc $IDAES_EXT/dist-petsc
fi
cp -r petscpy $IDAES_EXT/dist-petsc
cp ../dist-solvers/license.txt $IDAES_EXT/dist-petsc/license_petsc.txt
cp ../dist-solvers/version_solvers.txt $IDAES_EXT/dist-petsc/version_petsc.txt

# here you pack files
echo "#########################################################################"
echo "# Finish                                                                #"
echo "#########################################################################"
cd $IDAES_EXT/dist-petsc
tar -czvf idaes-petsc-${osname}-${MNAME}.tar.gz *
cd $IDAES_EXT/dist-solvers
tar -czvf idaes-solvers-${osname}-${MNAME}.tar.gz *
echo "Done"
echo "HSL Present: ${with_hsl}"
