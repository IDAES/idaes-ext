# Instructions

Currently you will need Windows to build the full set of binaries, since running
a Windows docker container requires Windows. The ```build-platforms``` directory
contains docker files to create build platform docker images.   The ```test-platforms```
directory contains docker files to create testing platform docker images. Each
directory contains a script ```build-windows.ps1``` to build a specific platform
run ```build-windows.ps1 <platform> --no-chache```.

Once you have test and build platform you can build and test the IDAES
extensions.  Docker files for build and testing are in the ```build-extensions```
and ```test-extensions``` directories.  Both directories contain an ```extras```
subdirectory, which contains files to be copied to the docker image.  For builds,
extras can contain HSL files.  For testing extras should contain the built binaries.
The builds and tests can be evoked with the ```build-windows.ps1 <platform> --no-chache```
script from within the appropriate directory.

The test containers install the appropriate binary package from extras with the
idaes-pse master branch and run tests.  The test containers can be run
interactively if you want to test against an IDAES repo that is not public.

To run a build or test container interactively do ```docker run -it ubuntu1910_test:latest```.
To delete the container created by the interactive session, list containers
with ```docker ps --all```. Then ```docker rm <CONTAINER ID>```.

To clean up orphaned images and containers ```docker system prune``` to delete an
image ```docker rmi <IMAGE>``` (e.g. ```docker rmi centos6_idaes_build```)
