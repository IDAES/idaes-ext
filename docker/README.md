# Docker build files

- `main-images/`: one Dockerfile per build flavor, plus `build-image.sh` to
  build the images.
- `build-scripts/`: `build-solvers.sh` (Linux) and `build.ps1` (Windows) build
  the tarballs inside those images.
  - `extras/`: put `coinhsl.tar.gz` here to build with HSL.
  - `tarballs/`: the built tarballs are written here.

See [doc/build.md](../doc/build.md) for the full build instructions.
