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
  
The powershell script will checkout a build container from DockerHub and execute the build scripts.  When the script completes, review the output for errors.  Source from outside repositories keeps changing, so it's not unusual to encounter errors. If there are errors, the solution is either to update the build scripts here or complain to the developers of the external pacakages. The binary tarball will be in the ```.\tarballs``` directory.  Repeat steps 3 and 6 for all platforms.

## Releasing a New Version

Once binaries have been compiled successfully for all platforms, they can be tested and released by the following procedure.  The eventual goal is to develop a GitHub Actions workflow to automatically do most of the release steps, but for now, many of the steps are manual.

1) Go to releases and locate the "Prerelease Binary Testing" release with the tag "test-release."
2) Edit the release to remove the tar-files from the previous release.
3) Add the newly created binary tar files to the "test-release" release.
4) Go to the actions tab, and manually run the Linux and Windows test actions.
5) Check the test results to ensure that the tests pass, or at least do not have any failures related to the binaries.
6) Delete the "new_release_testing" branch.
7) Run the "Release From Tests" action.  This will calculate hashes for the tar-files and commit a file "sha256_testing.txt" containing the hash information to the "releases" directory in the "new_release_testing" branch.
8) Rename the "sha256_testing.txt" by replacing testing with the tag to be used for the new release.
9) Create a PR from the "new_release_testing" branch to the "main" branch.  This should just add the new hash file, and nothing else.
10) Create the new release.  The tag should be the release number, for example "2.2.2".
11) Upload the binary tar-files as artifacts to the newly created release.  These should be the same ones tested previously, and should match the calculated hashes.
12) Update the default binary release tag in the IDAES-PSE repo.

