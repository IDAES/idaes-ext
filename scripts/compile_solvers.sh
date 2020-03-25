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
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
if [ -f $IDAES_EXT/../coinhsl.zip ]
then
  # If the HSL isn't there then just compile without.
  mkdir ThirdParty
  cd ThirdParty/
  # git clone -b stable/2.0 https://github.com/coin-or-tools/ThirdParty-HSL HSL
  mkdir HSL
  echo -n >HSL/.build
  #mkdir HSL/coinhsl
  #cp $IDAES_EXT/../coinhsl.zip HSL/coinhsl/
  #cd HSL/coinhsl
  #unzip coinhsl.zip
  cd $IDAES_EXT/coinbrew
else
  echo "HSL Not Available, BUILDING WITHOUT HSL" >&2
fi
bash coinbrew fetch $IPOPT_REPO:$IPOPT_BRANCH --no-prompt
bash coinbrew build $IPOPT_REPO:$IPOPT_BRANCH --no-prompt --disable-shared

cd $IDAES_EXT
mkdir dist-solvers
cd dist-solvers
cp ../coinbrew/dist/bin/ipopt ./
cp ../license.txt ./

# here you pack files
tar -czvf idaes-solvers-${osname}-64.tar.gz *
