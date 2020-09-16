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
bash coinbrew build Ipopt --no-prompt --disable-shared --enable-static LDFLAGS="-lgfortran -lm -llapack -lblas"

cd $IDAES_EXT
mkdir dist-solvers
cd dist-solvers
cp ../coinbrew/dist/bin/ipopt ./
cp ../license.txt ./
cp ../version.txt ./version_solvers.txt
sed s/"(DATE)"/`date +%Y%m%d-%H%M`/g version_solvers.txt > tmp
sed s/"(PLAT)"/${osname}/g tmp > tmp2
mv tmp2 version_solvers.txt


if [ "$(expr substr $(uname -s) 1 7)" == "MINGW64" ]
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

# Compile Pynumero

cd $IDAES_EXT
git clone $PYNU_REPO
patch pyomo/pyomo/contrib/pynumero/src/CMakeLists.txt < scripts/pynumero_cmake.diff
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
if [ "$(expr substr $(uname -s) 1 7)" == "MINGW64" ]
then
  cmake -DLINK_DLOPEN=OFF -DCMAKE_C_COMPILER=gcc -DCMAKE_LINKER=g++ -G"MSYS Makefiles" .
else
  cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_LINKER=g++ .
fi
make
cp bin/k_aug* $IDAES_EXT/dist-solvers

# here you pack files
cd $IDAES_EXT/dist-solvers
tar -czvf idaes-solvers-${osname}-64.tar.gz *
echo "HSL Present: ${with_hsl}"
