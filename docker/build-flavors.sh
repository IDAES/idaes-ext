#!/bin/bash

##
# Program usage
##

function usage() {
  printf "usage: $0 [flavor]\n"
  exit 0
}

##
# Choose flavor(s)
# Build all the flavors, or just the one specified
##

if [ $# -eq 0 ]; then
  dirs=$(find . -maxdepth 1 -type d -name "[a-zA-Z]*" -exec basename {} \;)
elif [ $# -eq 1 ]; then
  dirs="$1"
else
  usage
fi

##
# Build
##

# failure message
function ifailed() {
  msg=$1
  printf "*** FAILURE: $1 ***\n"
  exit 1
}

# status message
function imdoing() {
  printf "=== docker $1 ===\n"
}

# main loop
cid="x"
for d in ${dirs}; do
  printf "+=========================================\n"
  printf "|  $d\n"
  printf "+=========================================\n"
  cd $d || ifailed "no such directory: $d"
  # build
  imdoing build
  docker build --rm -t idaes:${d} . || ifailed "docker build"
  # run
  imdoing run
  cid=$(docker run -dt idaes:${d} /bin/bash)
  [ $? -eq 0 ] || ifailed "docker run"
  # cp
  imdoing cp
  cd ..
  docker cp ${cid}:/repo/idaes-ext/dist-lib/idaes-lib-${d}.tar.gz .
  if [ $? -ne 0 ]
  then
    imdoing stop
    docker stop ${cid}
    ifailed "docker cp"
  fi
  # done with this one
  imdoing stop
  docker stop ${cid}
  printf "=== SUCCESS: $d ===\n"
done

exit 0
