#!/bin/sh
# First argument is OS name
osname=$1; shift

# Make a directory to work in
export IDAES_EXT=`pwd`

# Set a few basic things

export IPOPT_BRANCH="idaes-3.13"
export IPOPT_REPO="https://github.com/idaes/Ipopt"
export PYNU_BRANCH="master"
export PYNU_REPO="https://github.com/pyomo/pyomo"
export K_AUG_BRANCH="ma57"
export K_AUG_REPO="https://github.com/dthierry/k_aug"

# Work-around for mumps gcc v10 gfortran bug
# Some reason this did not work on macOS. That is b/c gcc points to the clang/XCode version.
export GCC_VERSION=`gcc -dumpversion`
if [ "$(expr substr "$GCC_VERSION" 1 2)" = "10" ]; then
  export FCFLAGS="-w -fallow-argument-mismatch -O2"
  export FFLAGS="-w -fallow-argument-mismatch -O2"
fi

mkdir coinbrew
cd coinbrew
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew

# a line in the coinbrew script throws a syntax error on centos6 so change it
if [ ${osname} == "centos6" ]
then
  sed -e "s/\[ -v LD_LIBRARY_PATH \]/ \"1\" = \"0\" /g" -i coinbrew
fi

bash coinbrew fetch $IPOPT_REPO@$IPOPT_BRANCH --no-prompt
if [ -f $IDAES_EXT/../coinhsl.zip ]
then
  # if the HSL source zip is in place...
  echo -n >ThirdParty/HSL/.build
  
  # original
  mkdir ThirdParty/HSL/coinhsl
  cp $IDAES_EXT/../coinhsl.zip ThirdParty/HSL/coinhsl/
  cd ThirdParty/HSL/coinhsl
  
  # changes to get this working on macOS
  # cp $IDAES_EXT/../coinhsl.zip ThirdParty/HSL/
  # cd ThirdParty/HSL/
  unzip coinhsl.zip
  
  cd $IDAES_EXT/coinbrew
  with_hsl="YES"
else
  # If the HSL isn't there, build without it.
  echo "HSL Not Available, BUILDING SOLVERS WITHOUT HSL" >&2
  with_hsl="NO"
fi

# original
# bash coinbrew build Ipopt --no-prompt --disable-shared --enable-static LDFLAGS="-lgfortran -lm -llapack -lblas"

# adowling2 desktop
if [$(with_hsl) == "NO"]
then
  echo "Configuring coinbrew using --without-hsl flag"
  bash coinbrew build Ipopt --no-prompt --disable-shared --enable-static LDFLAGS="-lgfortran -lm -llapack -lblas -lgcc" --reconfigure CC="gcc-9" CXX="g++-9" F77="gfortran-9" --without-hsl
else
  echo "This should compile with HSL..."
  bash coinbrew build Ipopt --no-prompt --disable-shared --enable-static LDFLAGS="-lgfortran -lm -llapack -lblas -lgcc" --reconfigure CC="gcc-9" CXX="g++-9" F77="gfortran-9"
fi


# adowling2 desktop (with HSL)
#  bash coinbrew build Ipopt --no-prompt --disable-shared --enable-static LDFLAGS="-lgfortran -lm -llapack -lblas" --reconfigure CC="gcc-9" CXX="g++-9" F77="gfortran-9"

# adowling2 laptop
# bash coinbrew build Ipopt --no-prompt --disable-shared --enable-static LDFLAGS="-lgfortran -lm -llapack -lblas -lgcc" --reconfigure CC="gcc-10" CXX="g++-10" F77="gfortran-10" FCFLAGS="-w -fallow-argument-mismatch -O2" FFLAGS="-w -fallow-argument-mismatch -O2"

cd $IDAES_EXT
mkdir dist-solvers
cd dist-solvers
cp ../coinbrew/dist/bin/ipopt ./
cp ../license.txt ./
cp ../version.txt ./version_solvers.txt
sed s/"(DATE)"/`date +%Y%m%d-%H%M`/g version_solvers.txt > tmp
sed s/"(PLAT)"/${osname}/g tmp > tmp2
mv tmp2 version_solvers.txt


if [ "$(expr substr $(uname -s) 1 7)" = "MINGW64" ]
then
    # Winodws MinGW linked libraries
    cp /mingw64/bin/libstdc++-6.dll ./
    cp /mingw64/bin/libgcc_s_seh-1.dll ./
    cp /mingw64/bin/libwinpthread-1.dll ./
    cp /mingw64/bin/libgfortran-*.dll ./
    cp /mingw64/bin/libquadmath-0.dll ./
    cp /mingw64/bin/libgomp-1.dll ./
    cp /mingw64/bin/liblapack.dll ./
    cp /mingw64/bin/libblas.dll ./
fi

#if [ "$(expr substr $(uname -s) 1 6)" = "Darwin" ]
if [ "$(uname -s)" = "Darwin" ]
then
    echo "Copying libraries for macOS"
    # macOS linked libraries
    # For some reason, this code block is not executing...
    cp /usr/lib/libstdc++.6.dylib ./
    #cp /mingw64/bin/libgcc_s_seh-1.dll ./
    #cp /mingw64/bin/libwinpthread-1.dll ./
    #cp /mingw64/bin/libgfortran-*.dll ./
    cp /usr/local/opt/gcc/lib/gcc/9/libquadmath.0.dylib ./
    #cp /mingw64/bin/libgomp-1.dll ./
    #cp /mingw64/bin/liblapack.dll ./
    #cp /mingw64/bin/libblas.dll ./
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
  # Should this be changes to "gcc-9", etc.?
  cmake .. -DENABLE_HSL=no -DIPOPT_DIR=$IDAES_EXT/coinbrew/dist -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++
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

  # This is the original version
  # cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER.
  
  # This is my hack to get macOS to work
  # cmake -DCMAKE_C_COMPILER=gcc-9 .
  cmake -DCMAKE_C_COMPILER=gcc-10 .

fi
make
cp bin/k_aug* $IDAES_EXT/dist-solvers
cp dot_sens* $IDAES_EXT/dist-solvers

# here you pack files
cd $IDAES_EXT/dist-solvers
tar -czvf idaes-solvers-${osname}-64.tar.gz *
echo "HSL Present: ${with_hsl}"
