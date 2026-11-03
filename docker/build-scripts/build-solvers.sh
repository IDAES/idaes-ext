#!/usr/bin/env bash
# Build tarballs (written to tarballs/) in a Linux build image.
# Pass "" to use the default repo/branch/image. Extra args go to compile_solvers.sh.
# Put coinhsl.tar.gz in extras/ to build with HSL.
# Uses podman if available, else docker (override: CONTAINER_ENGINE).

set -euo pipefail

FLAVORS="el8 el9 ubuntu2004 ubuntu2204 ubuntu2404 ubuntu2604"

usage() {
  echo "Usage: $(basename "$0") <flavor> <arch> [repo] [branch] [image] [compile_solvers args...]"
  echo "  flavor: one of {${FLAVORS// /, }}"
  echo "  arch:   x86_64 or aarch64"
}

if [ $# -lt 2 ]; then
  usage
  exit 1
fi

flavor=$1
mname=$2
repo=${3:-https://github.com/IDAES/idaes-ext.git}
branch=${4:-main}
image=${5:-idaes-ext-${flavor}:${mname}}
if [ $# -gt 5 ]; then
  shift 5
  solver_buildargs="$*"
else
  solver_buildargs=""
fi

if [ "$flavor" = "windows" ]; then
  echo "windows is not supported; use build.ps1 on a Windows host." >&2
  exit 1
fi

case " $FLAVORS " in
  *" $flavor "*) ;;
  *) echo "Unknown flavor: $flavor" >&2; usage; exit 1 ;;
esac

case "$mname" in
  x86_64)  platform="linux/amd64" ;;
  aarch64) platform="linux/arm64" ;;
  *) echo "Unknown arch: $mname" >&2; usage; exit 1 ;;
esac

if [ -n "${CONTAINER_ENGINE:-}" ]; then
  engine=$CONTAINER_ENGINE
elif command -v podman >/dev/null 2>&1; then
  engine=podman
elif command -v docker >/dev/null 2>&1; then
  engine=docker
else
  echo "Neither podman nor docker was found." >&2
  exit 1
fi

wdir="/repo"
container="idaes-ext-${flavor}-${mname}-build-tmp"
here="$(cd "$(dirname "$0")" && pwd)"

echo "build-solvers.sh arguments:
    OS: $flavor
    Arch: $mname ($platform)
    Repo: $repo
    Branch: $branch
    Image: $image
    Engine: $engine
    Args to compile_solvers.sh: $solver_buildargs
"

cleanup() {
  "$engine" rm -f "$container" >/dev/null 2>&1 || true
}
trap cleanup EXIT
cleanup

"$engine" run --platform "$platform" --name "$container" -dt "$image"
"$engine" cp "$here/extras/" "$container:$wdir"
"$engine" exec "$container" sh -c "cp ${wdir}/extras/* ${wdir}"
"$engine" exec "$container" sh -c "cd ${wdir} && git clone ${repo} && cd idaes-ext && git checkout ${branch}"
"$engine" exec "$container" sh -c "cd ${wdir}/idaes-ext && bash scripts/compile_solvers.sh ${flavor} ${solver_buildargs}"
"$engine" exec "$container" sh -c "cd ${wdir}/idaes-ext && bash scripts/compile_libs.sh ${flavor}"
"$engine" stop "$container"

mkdir -p "$here/tarballs"
"$engine" cp "$container:$wdir/idaes-ext/dist-functions/idaes-functions-${flavor}-${mname}.tar.gz" "$here/tarballs/"
"$engine" cp "$container:$wdir/idaes-ext/dist/idaes-solvers-${flavor}-${mname}.tar.gz" "$here/tarballs/"

echo "Tarballs written to $here/tarballs/"
