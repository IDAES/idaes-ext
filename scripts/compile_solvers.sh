#!/bin/sh
# First argument is OS name
osname=$1; shift

# Make a directory to work in
export IDAES_EXT=`pwd`

# Set a few basic things

export IPOPT_BRANCH="idaes-3.13"
export IPOPT_REPO="https://github.com/idaes/Ipopt"

mkdir coinbrew
cd coinbrew
wget https://raw.githubusercontent.com/idaes/coinbrew/idaes-stable/coinbrew
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
  # If the HSL isn't there then just compile without.
  echo "HSL Not Available, NOT BUILDING SOLVERS" >&2
  exit 0
fi
bash coinbrew build Ipopt --no-prompt --disable-shared --enable-static LDFLAGS="-lgfortran -lm -llapack -lblas"

cd $IDAES_EXT
mkdir dist-solvers
cd dist-solvers
cp ../coinbrew/dist/bin/ipopt ./
cp ../license.txt ./

# here you pack files
tar -czvf idaes-solvers-${osname}-64.tar.gz *
