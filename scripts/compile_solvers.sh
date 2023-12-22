#!/bin/sh

# Exit on error.
set -e

# These scripts are not meant to be used for general builds.  They
# are taylor only to specific build systems.  They may provide some
# hints on how to build the solvers, but are limited.

# Set to exit on error
set -e

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
# Note that, when accessing patch files, we assume this directory
# to be the root of the idaes-ext repository.
export IDAES_EXT=`pwd`

arg2=$2
if [ $arg2 = "--without-hsl" ]; then
  echo "--without-hsl flag received. Building solvers without HSL."
  echo "HSL: NO"
  hslflag="--without-hsl"
  with_hsl="NO"
  build_hsl="NO"
elif [ -f $IDAES_EXT/../coinhsl.zip ]; then
  echo "coinhsl.zip found. Building solvers with HSL."
  echo "HSL: YES"
  hslflag="--with-hsl"
  with_hsl="YES"
  build_hsl="YES"
else
  echo "coinhsl.zip not found. Attempting to build with installed HSL."
  echo "HSL: YES"
  hslflag="--with-hsl"
  with_hsl="YES"
  build_hsl="NO"
fi

# Set a few basic things
export IPOPT_BRANCH="idaes-3.13"
export IPOPT_REPO="https://github.com/idaes/Ipopt"
export IPOPT_L1_BRANCH="restoration_mod"
export IPOPT_L1_REPO="https://github.com/idaes/Ipopt"
export PYNU_BRANCH="main"
export PYNU_REPO="https://github.com/pyomo/pyomo"
export K_AUG_BRANCH="default"
export K_AUG_REPO="https://github.com/idaes/k_aug"
export CC="gcc"
export CXX="g++"

# set PETSc location
# These petsc environment variables are used by the makefile
# in the $IDAES_EXT/petsc subdirectory.
if [ -z $PETSC_DIR ]; then
    # We only set PETSC_DIR if it is not already set. This is useful when
    # running this script locally (i.e. not via Docker)
    if [ ${osname} = "windows" ]
    then
      export PETSC_DIR=/c/repo/petsc-dist
    elif [ ${osname} = "darwin" ]; then
      export PETSC_DIR="$HOME/src/petsc-dist"
    else
      export PETSC_DIR=/repo/petsc-dist
    fi
    echo "PETSC_DIR has been set to $PETSC_DIR. If this is incorrect, set"
    echo "the PETSC_DIR environment variable before running this script."
fi
export PETSC_ARCH=""

# locate homebrew libs (this script is not at all for general builds)
if [ ${osname} = "darwin" ]
then
  if [ -f /opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib ]
  then
    export BREWLIB=/opt/homebrew/opt/gcc/lib/gcc/current/
  elif [ -f /usr/local/opt/gcc/lib/gcc/current/libgfortran.5.dylib ]; then
    export BREWLIB=/usr/local/opt/gcc/lib/gcc/current/
  fi
fi


mkdir coinbrew
cd coinbrew

if [ ${osname} = "darwin" ]; then
  curl --output coinbrew https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
  export CC="gcc-13"
  export CXX="g++-13"
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
  # We have already printed a message about HSL
  #echo "HSL is available, building with HSL"
  #with_hsl="YES"
#else
#  # If the HSL isn't there, build without it.
#  echo "HSL Not Available, BUILDING SOLVERS WITHOUT HSL" >&2
#  with_hsl="NO"
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
if [ $build_hsl = "YES" ]; then
  cd ThirdParty/HSL
  ./configure --disable-shared --enable-static --with-metis \
    --prefix=$IDAES_EXT/coinbrew/dist FFLAGS="-fPIC" CFLAGS="-fPIC" CXXFLAGS="-fPIC"
  make $PARALLEL
  make install
  cd $IDAES_EXT/coinbrew
fi

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
./configure --disable-shared --enable-static --with-mumps $hslflag \
  --prefix=$IDAES_EXT/coinbrew/dist
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

echo "#########################################################################"
echo "# Ipopt_L1 ampl executables                                             #"
echo "#########################################################################"
cd Ipopt_l1
if [ ${osname} = "el7" ]; then
  ./configure --disable-shared --enable-static --with-mumps $hslflag \
    ADD_CXXFLAGS="-std=c++11" \
    --prefix=$IDAES_EXT/coinbrew/dist_l1
else
  ./configure --disable-shared --enable-static --with-mumps $hslflag \
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
  --with-mumps $hslflag --enable-relocatable --prefix=$IDAES_EXT/coinbrew/dist-share
make $PARALLEL
make install
cd $IDAES_EXT/coinbrew

# The above compiles and installs solvers into ./coinbrew/dist*
# Now we will copy the files we wish to distribute into a ./dist directory,
# with standard bin, include, and lib sub-directories.

echo "#########################################################################"
echo "# Copy Coin Solver Files to dist                                        #"
echo "#########################################################################"
cd $IDAES_EXT
mkdir dist dist/bin dist/include dist/lib
# I'll try to avoid cd-ing around, for now, to make this script more
# explicit for the reader.
# cd dist
# Executables
if [ ${osname} = "windows" ]; then
  # windows
  cp ./coinbrew/dist_l1/bin/ipopt.exe ./dist/bin/ipopt_l1.exe
  cp ./coinbrew/dist_l1/bin/ipopt_sens.exe ./dist/bin/ipopt_sens_l1.exe
  # Explicitly only get ipopt so we don't get anything we shouldn't
  cp ./coinbrew/dist-share/bin/libipopt*.dll ./dist/lib/
  cp ./coinbrew/dist-share/bin/libsipopt*.dll ./dist/lib/
elif [ ${osname} = "darwin" ]; then
  cp ./coinbrew/dist_l1/bin/ipopt ./dist/bin/ipopt_l1
  cp ./coinbrew/dist_l1/bin/ipopt_sens ./dist/bin/ipopt_sens_l1
  # Explicitly only get ipopt so we don't get anything we shouldn't
  cp ./coinbrew/dist-share/lib/libipopt*.dylib ./dist/lib/
  cp ./coinbrew/dist-share/lib/libsipopt*.dylib ./dist/lib/
else
  # linux
  cp ./coinbrew/dist_l1/bin/ipopt ./dist/bin/ipopt_l1
  cp ./coinbrew/dist_l1/bin/ipopt_sens ./dist/bin/ipopt_sens_l1
  # Explicitly only get ipopt so we don't get anything we shouldn't
  cp ./coinbrew/dist-share/lib/libipopt*.so ./dist/lib/
  cp ./coinbrew/dist-share/lib/libsipopt*.so ./dist/lib/
fi
cp ./coinbrew/dist/bin/ipopt ./dist/bin/
cp ./coinbrew/dist/bin/ipopt_sens ./dist/bin/
cp ./coinbrew/dist/bin/clp ./dist/bin/
cp ./coinbrew/dist/bin/cbc ./dist/bin/
cp ./coinbrew/dist/bin/bonmin ./dist/bin/
cp ./coinbrew/dist/bin/couenne ./dist/bin/
# Run strip to remove unneeded symbols (it's okay that some files
#  aren't exe or libraries)
if [ $osname = "darwin" ]; then
  # I'm not sure if this exactly reproduces --strip-unneeded on darwin
  strip ./dist/bin/*
else
  strip --strip-unneeded ./dist/bin/*
fi

# Copy contents of dist-share (ipopt lib, include, and share) into dist
cp -r ./coinbrew/dist-share/* ./dist/

# Patch to remove "Requires.private: coinhsl coinmumps"
# Note that we are patching the file *after* we copy it into the directory
# we will distribute.
#patch ./dist/lib/pkgconfig/ipopt.pc < $IDAES_EXT/scripts/ipopt.pc.patch

# Handle platform-dependent sed behavior...
if [ $osname = "darwin" ]; then
  sed -i "" "s/^Requires\.private/#Requires.private/" dist/lib/pkgconfig/ipopt.pc
else
  sed -i "s/^Requires\.private/#Requires.private/" dist/lib/pkgconfig/ipopt.pc
fi

echo "#########################################################################"
echo "# Copy License and Version Files to dist-solvers                        #"
echo "#########################################################################"
# Text information files include build time
cp ./license.txt ./dist/
cp ./version.txt ./dist/version_solvers.txt
sed s/"(DATE)"/`date +%Y%m%d-%H%M`/g ./dist/version_solvers.txt > ./dist/tmp
sed s/"(PLAT)"/${osname}-${MNAME}/g ./dist/tmp > ./dist/tmp2
# Why do we create dist/version_solvers.txt above if we will just overwrite it?
mv ./dist/tmp2 ./dist/version_solvers.txt
rm ./dist/tmp*

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
    cp /mingw64/bin/libstdc++-6.dll ./dist/lib/
    cp /mingw64/bin/libgcc_s_seh-1.dll ./dist/lib/
    cp /mingw64/bin/libwinpthread-1.dll ./dist/lib/
    cp /mingw64/bin/libgfortran-*.dll ./dist/lib/
    cp /mingw64/bin/libquadmath-0.dll ./dist/lib/
    cp /mingw64/bin/libgomp-1.dll ./dist/lib/
    cp /mingw64/bin/liblapack.dll ./dist/lib/
    cp /mingw64/bin/libblas.dll ./dist/lib/
    cp /mingw64/bin/libbz2-*.dll ./dist/lib/
    cp /mingw64/bin/zlib*.dll ./dist/lib/
    cp /mingw64/bin/libssp*.dll ./dist/lib/
fi

if [ ${osname} = "darwin" ]; then
  cp ${BREWLIB}libgfortran.5.dylib ./dist/lib/
  cp ${BREWLIB}libgcc_s.1.1.dylib ./dist/lib/
  cp ${BREWLIB}libstdc++.6.dylib ./dist/lib/
  cp ${BREWLIB}libgomp.1.dylib ./dist/lib/
  cp ${BREWLIB}libquadmath.0.dylib ./dist/lib/
fi

echo "#########################################################################"
echo "# Pynumero                                                              #"
echo "#########################################################################"

if [ ${osname} = "darwin" ]; then
  export CC="cc"
  export CXX="c++"
fi

# We are already in this directory
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
cp libpynumero_ASL* $IDAES_EXT/dist/lib/
# Return to root directory
cd $IDAES_EXT

echo "#########################################################################"
echo "# k_aug, dotsens                                                        #"
echo "#########################################################################"
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
cp bin/k_aug* $IDAES_EXT/dist/bin/
cp dot_sens* $IDAES_EXT/dist/bin/
# Return to root directory
cd $IDAES_EXT

echo "#########################################################################"
echo "# PETSc                                                                 #"
echo "#########################################################################"
export ASL_INC=$IDAES_EXT/coinbrew/dist/include/coin-or/asl
export ASL_LIB=$IDAES_EXT/coinbrew/dist/lib/libcoinasl.a
cd $IDAES_EXT/petsc
make $PARALLEL
make py
if [ ${osname} = "windows" ]
then
  cp petsc.exe $IDAES_EXT/dist/bin/
else
  cp petsc $IDAES_EXT/dist/bin/
fi
cp -r petscpy $IDAES_EXT/dist/lib/

# Return to root directory
cd $IDAES_EXT

if [ ${osname} = "darwin" ]; then
  echo "#########################################################################"
  echo "# macOS update rpaths                                                   #"
  echo "#########################################################################"
  update_rpath_darwin() {
    install_name_tool -change ${BREWLIB}libgfortran.5.dylib @rpath/libgfortran.5.dylib $1
    install_name_tool -change ${BREWLIB}libgcc_s.1.1.dylib @rpath/libgcc_s.1.1.dylib $1
    install_name_tool -change ${BREWLIB}libstdc++.6.dylib @rpath/libstdc++.6.dylib $1
    install_name_tool -change ${BREWLIB}libgomp.1.dylib @rpath/libgomp.1.dylib $1
    install_name_tool -change ${BREWLIB}libquadmath.0.dylib @rpath/libquadmath.0.dylib $1
  }
  # Dynamic libraries also need to update their own install names
  update_library_rpath_darwin() {
    install_name_tool -id @rpath/$1 $1
  }
  update_rpath_darwin dist/bin/ipopt
  update_rpath_darwin dist/bin/ipopt_sens
  update_rpath_darwin dist/bin/clp
  update_rpath_darwin dist/bin/cbc
  update_rpath_darwin dist/bin/bonmin
  update_rpath_darwin dist/bin/couenne  
  update_rpath_darwin dist/bin/ipopt_l1
  update_rpath_darwin dist/bin/ipopt_sens_l1
  update_rpath_darwin dist/lib/libipopt.dylib
  update_rpath_darwin dist/lib/libsipopt.dylib
  update_rpath_darwin dist/lib/libipopt.3.dylib
  update_rpath_darwin dist/lib/libsipopt.3.dylib
  update_rpath_darwin dist/lib/libpynumero_ASL.dylib
  update_library_rpath_darwin dist/lib/libipopt.dylib
  update_library_rpath_darwin dist/lib/libsipopt.dylib
  update_library_rpath_darwin dist/lib/libipopt.3.dylib
  update_library_rpath_darwin dist/lib/libsipopt.3.dylib
  update_library_rpath_darwin dist/lib/libpynumero_ASL.dylib
  # if no hsl k_aug and dot_snse won't exist
  update_rpath_darwin dist/bin/k_aug || true 
  update_rpath_darwin dist/bin/dot_sens || true
  update_rpath_darwin dist/bin/petsc
fi

# here you pack files
echo "#########################################################################"
echo "# Finish                                                                #"
echo "#########################################################################"
cd dist
tar -czvf idaes-solvers-${osname}-${MNAME}.tar.gz *
cd $IDAES_EXT
echo "Done"
echo "HSL Present: ${with_hsl}"
