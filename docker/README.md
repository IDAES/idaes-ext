# Instructions

This directory can build Docker containers that compile the
binary extensions for various Linux flavors.

To build all flavors (sequentially)

    # ./build-flavors.sh

To build one flavor directory, e.g. "centos7"

    # ./build-flavors.sh centos7

Look for "=== SUCCESS" or "*** FAILURE" at the bottom of the output.

If it succeeeds, the result(s) will be files in this
directory, named for each flavor. For example, for centos7 you will
get `idaes-lib-centos7.tar.gz`.
