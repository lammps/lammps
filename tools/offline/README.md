# Building LAMMPS and its documentation on offline systems

In some situations it might be necessary to build LAMMPS on a system without
internet. The scripts in this folder allow you to preload external dependencies
for both the documentation build for building with CMake into a folder and then
use that folder on an offline system.

It does so by

1.) Downloading necessary pip packages
2.) Cloning git repositories
3.) Downloading tarballs

As of April 2021, all of these downloads make up around 600MB.  By
default, it will download everything into $HOME/.cache/lammps, but this can be
changed with the ``LAMMPS_CACHING_DIR`` environment variable.

Once the caches have been initialized, they can be used for building
LAMMPS documentation or compiling using CMake on an offline system.

The ``use_caches.sh`` must be sourced into the current shell to initialize the
offline build environment. Note that it must use the same ``LAMMPS_CACHING_DIR``.
This script does the following:

1.) Sets up environment variables that modify the behavior of both pip and git
2.) Starts a simple local HTTP server to host files for CMake

Afterwards, it will print out instruction on how to modify the CMake command
line to make sure it uses the local HTTP server.

To undo the environment changes and shutdown the HTTP server, run the
``deactivate_caches`` command.

## Examples

For all of the examples below, you first need to create the cache (which requires internet).

```bash
./tools/offline/init_caches.sh
```

Afterwards, you can disconnect or copy the contents of the
``LAMMPS_CACHING_DIR`` folder to an offline system.

### Documentation

```bash
# if LAMMPS_CACHING_DIR is different from default, make sure to set it first
# export LAMMPS_CACHING_DIR=path/to/folder
source tools/offline/use_caches.sh
cd doc/
make html

deactivate_caches
```

### CMake Build

```bash
# if LAMMPS_CACHING_DIR is different from default, make sure to set it first
# export LAMMPS_CACHING_DIR=path/to/folder
source tools/offline/use_caches.sh

mkdir build
cd build
cmake -D LAMMPS_DOWNLOADS_URL=${HTTP_CACHE_URL} -C ${LAMMPS_HTTP_CACHE_CONFIG} -C ../cmake/presets/most.cmake ../cmake
make -j 8

deactivate_caches
```
