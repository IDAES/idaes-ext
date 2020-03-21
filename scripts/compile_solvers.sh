#!/bin/sh
# First argument is OS name
osname=$1; shift

# Make a directory to work in
export IDAES_EXT=`pwd`

# Set a few basic things

export IPOPT_BRANCH="idaes-3.13.0"
export IPOPT_REPO="https://github.com/idaes/Ipopt"

mkdir coinbrew
cd coinbrew
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
bash coinbrew fetch $IPOPT_REPO:$IPOPT_BRANCH --no-prompt
if [ -d $IDAES_EXT/coinhsl ]
then
  # If the HSL isn't there then just compile without.
  echo -n >ThirdParty/HSL/.build
  cp -r $IDAES_EXT/coinhsl ThirdParty/HSL/coinhsl
else
  echo "HSL Not Available, BUILDING WITHOUT HSL" >&2
fi
bash coinbrew build $IPOPT_REPO:$IPOPT_BRANCH --no-prompt

cd $IDAES_EXT
mkdir dist-solvers
cd dist-solvers
cp ../coinbrew/dist/bin/ipopt ./
cp ../license.txt ./

# here you pack files
tar -czvf idaes-solvers-${osname}-64.tar.gz *
