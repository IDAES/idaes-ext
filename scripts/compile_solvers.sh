#!/bin/sh
# First argument is OS name
osname=$1; shift

# Make a directory to work in
export IDAES_EXT=`pwd`

# Set a few basic things

export IPOPT_BRANCH="idaes-3.13"
export IPOPT_REPO="https://github.com/idaes/Ipopt"
export PYNU_BRANCH="pynumero-testing"
export PYNU_REPO="https://github.com/jsiirola/pyomo"

mkdir coinbrew
cd coinbrew
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew

# a line in the coinbrew script throws a syntax error on centos6 so change it
if [ ${osname}=="centos6" ]
then
  sed -e "s/\[ -v LD_LIBRARY_PATH \]/ \"1\" = \"0\" /g" -i coinbrew
fi

bash coinbrew fetch $IPOPT_REPO:$IPOPT_BRANCH --no-prompt
if [ -f $IDAES_EXT/../coinhsl.zip ]
then
  # if the HSL source zip is in place...
  echo -n >ThirdParty/HSL/.build
  mkdir ThirdParty/HSL/coinhsl
  cp $IDAES_EXT/../coinhsl.zip ThirdParty/HSL/coinhsl/
  cd ThirdParty/HSL/coinhsl
  unzip coinhsl.zip
  cd $IDAES_EXT/coinbrew
else
  # If the HSL isn't there then just quit.  We could also build without if
  # there is a reason.
  echo "HSL Not Available, NOT BUILDING SOLVERS" >&2
  exit 0
fi
bash coinbrew build Ipopt --no-prompt --disable-shared --enable-static LDFLAGS="-lgfortran -lm -llapack -lblas"

cd $IDAES_EXT
mkdir dist-solvers
cd dist-solvers
cp ../coinbrew/dist/bin/ipopt ./
cp ../license.txt ./
cp ../version.txt ./
cp ../amplmp_license.rst ./

if [ "$(expr substr $(uname -s) 1 7)" == "MINGW64" ]
then
    # Winodws MinGW linked libraries
    cp /mingw64/bin/libstdc++-6.dll ./
    cp /mingw64/bin/libgcc_s_seh-1.dll ./
    cp /mingw64/bin/libwinpthread-1.dll ./
    cp /mingw64/bin/libgfortran-4.dll ./
    cp /mingw64/bin/libquadmath-0.dll ./
    cp /mingw64/bin/libgomp-1.dll ./
    cp /mingw64/bin/liblapack.dll ./
    cp /mingw64/bin/libblas.dll ./
fi

# Compile Pynumero

cd $IDAES_EXT
git clone $PYNU_REPO
cd pyomo
git checkout $PYNU_BRANCH
cd pyomo/contrib/pynumero/src
mkdir build
cd build
if [ "$(expr substr $(uname -s) 1 7)" == "MINGW64" ]
then
  cmake -G"MSYS Makefiles" .. -DBUILD_AMPLMP_IF_NEEDED=ON
else
  cmake .. -DBUILD_AMPLMP_IF_NEEDED=ON
fi
make
cp libpynumero_ASL* $IDAES_EXT/dist-solvers
#cp sparse_utils/libpynumero_SPARSE* $IDAES_EXT/dist-lib
cd $IDAES_EXT/dist-solvers

# here you pack files
tar -czvf idaes-solvers-${osname}-64.tar.gz *
