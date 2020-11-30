# IDAES Extensions

This repository hosts IDAES binary extensions.  These extensions include solvers and function libraries.  These binary extensions are intended for use with the IDAES Process Modeling Framework https://github.com/IDAES/idaes-pse.

## Testing

Latest binary test results are [here](test_status.md). The results are just for the latest testing release.  All othere released binaries were tested and work with the corresonding IDEAS-PSE release.

## Getting Extensions

If you have the IDAES framework installed, you can get the extensions by running the following command:

```sh
idaes get-extensions
```

You can get a particular release tag (2.1.1 and higher) by providing a release tag:

```sh
idaes get-extensions --release <tag>
```

You can also download a particular release of the extensions by providing a URL:

```sh
idaes get-extensions --url <url>
```

To get a build for a particular platform, there is a platform option. The get-extensions command can detect operating system and common Linux distributions, so specifying a platform is rarely required. Current platform options are: centos6, centos7, centos8, rhel6, rhel7, rhel8, ubuntu1804, ubuntu1910, ubuntu2004, windows, and darwin. Darwin (Mac) is not currently available, but will be soon.  Currently, if no platform is specified or detected on Linux, the "centos7" build will be used. If you are using an unsuported Linux distribution, specifying a similar platform will likely work anyway.

```sh
idaes get-extensions --platform <platform>
```

## Contents

The extensions contain a version of the IPOPT solver compiled the HSL linear solver library, for which IDAES has obtained
a distribution license. All technical papers, sales and publicity material resulting from use of the HSL codes within IPOPT
must contain the following acknowledgement: HSL, a collection of Fortran codes for large-scale scientific computation. See http://www.hsl.rl.ac.uk. The Ipopt source code for IDAES was modified only to add additional messages as required by the HSL
distribution license.  The modified Ipopt code is available here https://github.com/IDAES/Ipopt/tree/idaes-3.12.13.

A build Pynumero-ASL is included, see https://github.com/Pyomo/pyomo/tree/master/pyomo/contrib/pynumero.  The Pynumero build specifically excludes the HSL functions.

A build of k\_aug is included that used the HSL ma57 linear solver, for which the IDAES project has obtained a distribution license.  For more information on, k\_aug see https://github.com/dthierry/k_aug.   

See https://github.com/IDAES/idaes-ext/blob/master/license.txt for additional information on third-party code contained
in the binaries.

## Building

These build instructions require Windows since Windows can run both Linux and Windows docker containers.  You will also need to have Docker Desktop installed.

1) Download the idaes-ext repository as a zip-file.  We suggest not cloning the repo when doing a build to avoid inadvertantly pushing back changes or files that shouldn't be there.
2) Unzip the idaes-ext zip-file in a convenient location.
3) Make sure Docker Desktop is in the appropriate mode, Linux or Windows, for the binaries to be built.
4) In PowerShell, change to the directory ```idaes-ext\docker\build-extensions```
5) Copy ```coinhsl.zip``` into the ```.\extras``` directory.
6) Execute the build script, where *platform* is in {windows, ubuntu1804, ubuntu1910, ubuntu2004, centos6, centos7, centos8}.
  > .\build-windows.ps1 *platform* --no-cache
  
The build scripts will checkout a build container from DockerHub and execute the build scripts.  When the script completes, review the output for errors.  Source from outside repositories keeps changing, so it's not unusual to encounter errors. If there are errors, the solution is either to update the build scripts here or complain to the developers of the external pacakages. The binary tarball will be in the ```.\tarballs``` directory.  Repeat steps 3 and 6 for all platforms.

## Releasing a New Version

Todo
