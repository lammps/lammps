![Kokkos](https://avatars2.githubusercontent.com/u/10199860?s=200&v=4)

# Kokkos: Core Libraries

Kokkos Core implements a programming model in C++ for writing performance portable
applications targeting all major HPC platforms. For that purpose it provides
abstractions for both parallel execution of code and data management.
Kokkos is designed to target complex node architectures with N-level memory
hierarchies and multiple types of execution resources. It currently can use
CUDA, HPX, OpenMP and Pthreads as backend programming models with several other
backends in development.

Kokkos Core is part of the Kokkos C++ Performance Portability Programming EcoSystem,
which also provides math kernels (https://github.com/kokkos/kokkos-kernels), as well as 
profiling and debugging tools (https://github.com/kokkos/kokkos-tools).  

# Learning about Kokkos

A programming guide can be found on the Wiki, the API reference is under development.

For questions find us on Slack: https://kokkosteam.slack.com or open a github issue.

For non-public questions send an email to
crtrott(at)sandia.gov

A separate repository with extensive tutorial material can be found under 
https://github.com/kokkos/kokkos-tutorials.

Furthermore, the 'example/tutorial' directory provides step by step tutorial
examples which explain many of the features of Kokkos. They work with
simple Makefiles. To build with g++ and OpenMP simply type 'make'
in the 'example/tutorial' directory. This will build all examples in the
subfolders. To change the build options refer to the Programming Guide
in the compilation section.

To learn more about Kokkos consider watching one of our presentations:
* GTC 2015:
  - http://on-demand.gputechconf.com/gtc/2015/video/S5166.html
  - http://on-demand.gputechconf.com/gtc/2015/presentation/S5166-H-Carter-Edwards.pdf


# Contributing to Kokkos

We are open and try to encourage contributions from external developers. 
To do so please first open an issue describing the contribution and then issue
a pull request against the develop branch. For larger features it may be good
to get guidance from the core development team first through the github issue. 

Note that Kokkos Core is licensed under standard 3-clause BSD terms of use. 
Which means contributing to Kokkos allows anyone else to use your contributions
not just for public purposes but also for closed source commercial projects.
For specifics see the LICENSE file contained in the repository or distribution.

# Requirements

### Primary tested compilers on X86 are:
* GCC 4.8.4
* GCC 4.9.3
* GCC 5.1.0
* GCC 5.4.0
* GCC 5.5.0
* GCC 6.1.0
* GCC 7.2.0
* GCC 7.3.0
* GCC 8.1.0
* Intel 15.0.2
* Intel 16.0.1
* Intel 17.0.1
* Intel 17.4.196
* Intel 18.2.128
* Clang 3.6.1
* Clang 3.7.1
* Clang 3.8.1
* Clang 3.9.0
* Clang 4.0.0
* Clang 6.0.0 for CUDA (CUDA Toolkit 9.0)
* Clang 7.0.0 for CUDA (CUDA Toolkit 9.1)
* Clang 8.0.0 for CUDA (CUDA Toolkit 9.2)
* PGI 18.7
* NVCC 9.1 for CUDA (with gcc 6.1.0)
* NVCC 9.2 for CUDA (with gcc 7.2.0)
* NVCC 10.0 for CUDA (with gcc 7.4.0)
* NVCC 10.1 for CUDA (with gcc 7.4.0)

### Primary tested compilers on Power 8 are:
* GCC 6.4.0 (OpenMP,Serial)
* GCC 7.2.0 (OpenMP,Serial)
* IBM XL 16.1.0 (OpenMP, Serial)
* NVCC 9.2.88 for CUDA (with gcc 7.2.0 and XL 16.1.0)

### Primary tested compilers on Intel KNL are:
* Intel 16.4.258 (with gcc 4.7.2)
* Intel 17.2.174 (with gcc 4.9.3)
* Intel 18.2.199 (with gcc 4.9.3)

### Primary tested compilers on ARM (Cavium ThunderX2)
* GCC 7.2.0 
* ARM/Clang 18.4.0
  
### Other compilers working:
* X86:
    * Cygwin 2.1.0 64bit with gcc 4.9.3
    * GCC 8.1.0 (not warning free)

### Known non-working combinations:
* Power8:
    * Pthreads backend
* ARM
    * Pthreads backend


Primary tested compiler are passing in release mode
with warnings as errors. They also are tested with a comprehensive set of 
backend combinations (i.e. OpenMP, Pthreads, Serial, OpenMP+Serial, ...).
We are using the following set of flags:
* GCC:   
   ````
      -Wall -Wshadow -pedantic 
      -Werror -Wsign-compare -Wtype-limits
      -Wignored-qualifiers -Wempty-body 
      -Wclobbered -Wuninitialized
   ````
* Intel: 
    ````
      -Wall -Wshadow -pedantic 
      -Werror -Wsign-compare -Wtype-limits 
      -Wuninitialized
    ````
* Clang: 
    ````
      -Wall -Wshadow -pedantic 
      -Werror -Wsign-compare -Wtype-limits 
      -Wuninitialized
    ````    

* NVCC:  
  ````
    -Wall -Wshadow -pedantic 
    -Werror -Wsign-compare -Wtype-limits 
    -Wuninitialized
  ````

Other compilers are tested occasionally, in particular when pushing from develop to 
master branch. These are tested less rigorously without `-Werror` and only for a select set of backends.

# Building and Installing Kokkos
Kokkos provide a CMake build system and a raw Makefile build system. 
The CMake build system is strongly encouraged and will be the most rigorously supported in future releases.
Full details are given in the [build instructions](BUILD.md). Basic setups are shown here:

## CMake

The best way to install Kokkos is using the CMake build system. Assuming Kokkos lives in `$srcdir`: 
````
cmake $srcdir \
  -DCMAKE_CXX_COMPILER=$path_to_compiler \
  -DCMAKE_INSTALL_PREFIX=$path_to_install \
  -DKokkos_ENABLE_OPENMP=On \
  -DKokkos_ARCH_HSW=On \
  -DKokkos_ENABLE_HWLOC=On \
  -DKokkos_HWLOC_DIR=$path_to_hwloc
````
then simply type `make install`. The Kokkos CMake package will then be installed in `$path_to_install` to be used by downstream packages.

To validate the Kokkos build, configure with 
````
 -DKokkos_ENABLE_TESTS=On 
````
and run `make test` after completing the build.

For your CMake project using Kokkos, code such as the following:

````
find_package(Kokkos)
...
target_link_libraries(myTarget Kokkos::kokkos)
````
should be added to your CMakeLists.txt. Your configure should additionally include
````
-DKokkos_DIR=$path_to_install/cmake/lib/Kokkos
````
or
````
-DKokkos_ROOT=$path_to_install
````
for the install location given above.

## Spack
An alternative to manually building with the CMake is to use the Spack package manager.
To do so, download the `kokkos-spack` git repo and add to the package list:
````
spack repo add $path-to-kokkos-spack
````
A basic installation would be done as:
````
spack install kokkos
````
Spack allows options and and compilers to be tuned in the install command.
````
spack install kokkos@3.0 %gcc@7.3.0 +openmp
````
This example illustrates the three most common parameters to Spack:
* Variants: specified with, e.g. `+openmp`, this activates (or deactivates with, e.g. `~openmp`) certain options.
* Version:  immediately following `kokkos` the `@version` can specify a particular Kokkos to build
* Compiler: a default compiler will be chosen if not specified, but an exact compiler version can be given with the `%`option.

For a complete list of Kokkos options, run:
````
spack info kokkos
````
Spack currently installs packages to a location determined by a unique hash. This hash name is not really "human readable".
Generally, Spack usage should never really require you to reference the computer-generated unique install folder. 
More details are given in the [build instructions](BUILD.md). If you must know, you can locate Spack Kokkos installations with:
````
spack find -p kokkos ...
````
where `...` is the unique spec identifying the particular Kokkos configuration and version.


## Raw Makefile 
A bash script is provided to generate raw makefiles.
To install Kokkos as a library create a build directory and run the following
````
$KOKKOS_PATH/generate_makefile.bash --prefix=$path_to_install
````
Once the Makefile is generated, run:
````
make kokkoslib
make install
````
To additionally run the unit tests:
````
make build-test
make test
````
Run `generate_makefile.bash --help` for more detailed options such as
changing the device type for which to build.

## Inline Builds vs. Installed Package
For individual projects, it may be preferable to build Kokkos inline rather than link to an installed package.
The main reason is that you may otherwise need many different
configurations of Kokkos installed depending on the required compile time
features an application needs. For example there is only one default 
execution space, which means you need different installations to have OpenMP
or Pthreads as the default space. Also for the CUDA backend there are certain
choices, such as allowing relocatable device code, which must be made at 
installation time. Building Kokkos inline uses largely the same process
as compiling an application against an installed Kokkos library. 

For CMake, this means copying over the Kokkos source code into your project and adding `add_subdirectory(kokkos)` to your CMakeLists.txt.

For raw Makefiles, see the example benchmarks/bytes_and_flops/Makefile which can be used with an installed library and or an inline build.  

# Kokkos and CUDA UVM

Kokkos does support UVM as a specific memory space called CudaUVMSpace. 
Allocations made with that space are accessible from host and device. 
You can tell Kokkos to use that as the default space for Cuda allocations.
In either case UVM comes with a number of restrictions:
* You can't access allocations on the host while a kernel is potentially 
running. This will lead to segfaults. To avoid that you either need to 
call Kokkos::Cuda::fence() (or just Kokkos::fence()), after kernels, or
you can set the environment variable CUDA_LAUNCH_BLOCKING=1.
* In multi socket multi GPU machines without NVLINK, UVM defaults 
to using zero copy allocations for technical reasons related to using multiple
GPUs from the same process. If an executable doesn't do that (e.g. each
MPI rank of an application uses a single GPU [can be the same GPU for 
multiple MPI ranks]) you can set CUDA_MANAGED_FORCE_DEVICE_ALLOC=1.
This will enforce proper UVM allocations, but can lead to errors if 
more than a single GPU is used by a single process.


# Citing Kokkos

If you publish work which mentions Kokkos, please cite the following paper:

````
@article{CarterEdwards20143202,
  title = "Kokkos: Enabling manycore performance portability through polymorphic memory access patterns ",
  journal = "Journal of Parallel and Distributed Computing ",
  volume = "74",
  number = "12",
  pages = "3202 - 3216",
  year = "2014",
  note = "Domain-Specific Languages and High-Level Frameworks for High-Performance Computing ",
  issn = "0743-7315",
  doi = "https://doi.org/10.1016/j.jpdc.2014.07.003",
  url = "http://www.sciencedirect.com/science/article/pii/S0743731514001257",
  author = "H. Carter Edwards and Christian R. Trott and Daniel Sunderland"
}
````

##### [LICENSE](https://github.com/kokkos/kokkos/blob/master/LICENSE)

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Under the terms of Contract DE-NA0003525 with NTESS,
the U.S. Government retains certain rights in this software.

