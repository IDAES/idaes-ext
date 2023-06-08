# README

This repository contains scripts and source code for compiling IDAES extensions
these consist of IDAES property libraries and third party solvers. IDAES users
generally do not need to access this repo.  The build scripts are specific
to the IDAES build environment and not intended for general use.

More detailed build docs are here https://github.com/IDAES/idaes-ext/blob/main/doc/build.md.

For more information on building binaries see ./doc/build.md.  For more information
about the content of binary releases see https://idaes-pse.readthedocs.io/en/stable/tutorials/getting_started/binaries.html#binary-packages.
For installation instructions see https://idaes-pse.readthedocs.io/en/stable/tutorials/getting_started/index.html.

## Testing status

Currently, since the builds must be done off-line and uploaded to a release, the
test status does not necessarily reflect the status of the code currently in the
repository or any third-party repository.  This just shows the status of the most
recently tested builds.

### Windows

[![Windows Test IDAES-PSE Main](https://github.com/IDAES/idaes-ext/actions/workflows/test_windows_main.yml/badge.svg)](https://github.com/IDAES/idaes-ext/actions/workflows/test_windows_main.yml)

### Linux

Currently Testing:
* CentOS 7
* Rocky Linlux 8
* Ubuntu 18.04
* Ubuntu 20.04
* Ubuntu 22.04
* Debian 9
* Debian 10
* Debian 11

[![Linux Test IDAES-PSE Main](https://github.com/IDAES/idaes-ext/actions/workflows/test_linux.yml/badge.svg)](https://github.com/IDAES/idaes-ext/actions/workflows/test_linux.yml)

### macOS

#### Intel

[![macOS Test IDAES-PSE Main](https://github.com/IDAES/idaes-ext/actions/workflows/test_macos_main.yml/badge.svg)](https://github.com/IDAES/idaes-ext/actions/workflows/test_macos_main.yml)

#### Apple Silicon

Coming Soon.
