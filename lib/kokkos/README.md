![Kokkos](https://avatars2.githubusercontent.com/u/10199860?s=200&v=4)

# Kokkos: Core Libraries

Kokkos Core implements a programming model in C++ for writing performance portable
applications targeting all major HPC platforms. For that purpose it provides
abstractions for both parallel execution of code and data management.
Kokkos is designed to target complex node architectures with N-level memory
hierarchies and multiple types of execution resources. It currently can use
CUDA, HIP, SYCL, HPX, OpenMP and C++ threads as backend programming models with several other
backends in development.

Kokkos Core is part of the Kokkos C++ Performance Portability Programming EcoSystem,
which also provides math kernels (https://github.com/kokkos/kokkos-kernels), as well as
profiling and debugging tools (https://github.com/kokkos/kokkos-tools).

# Learning about Kokkos

The best way to start learning about Kokkos is going through the Kokkos Lectures.
They are online available at https://kokkos.link/the-lectures and contain a mix
of lecture videos and hands-on exercises covering all the important Kokkos Ecosystem
capabilities.

A programming guide and API reference can be found on the Wiki
(https://github.com/kokkos/kokkos/wiki).

For questions find us on Slack: https://kokkosteam.slack.com or open a github issue.

For non-public questions send an email to
crtrott(at)sandia.gov

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

### Minimum Compiler Versions

Generally Kokkos should work with all compiler versions newer than the minimum.
However as in all sufficiently complex enough code, we have to work around compiler
bugs with almost all compilers. So compiler versions we don't test may have issues
we are unaware of.

* GCC: 5.3.0
* Clang: 4.0.0
* Intel: 17.0.1
* NVCC: 9.2.88
* NVC++: 21.5
* ROCm: 4.3
* MSVC: 19.29
* IBM XL: 16.1.1
* Fujitsu: 4.5.0
* ARM/Clang 20.1

### Primary Tested Compilers

* GCC: 5.3.0, 6.1.0, 7.3.0, 8.3, 9.2, 10.0
* NVCC: 9.2.88, 10.1, 11.0
* Clang: 8.0.0, 9.0.0, 10.0.0, 12.0.0
* Intel 17.4, 18.1, 19.5
* MSVC: 19.29
* ARM/Clang: 20.1
* IBM XL: 16.1.1
* ROCm: 4.3.0

### Build system:

* CMake >= 3.16: required
* CMake >= 3.18: Fortran linkage. This does not affect most mixed Fortran/Kokkos builds. See [build issues](BUILD.md#KnownIssues).
* CMake >= 3.21.1 for NVC++

Primary tested compiler are passing in release mode
with warnings as errors. They also are tested with a comprehensive set of
backend combinations (i.e. OpenMP, Threads, Serial, OpenMP+Serial, ...).
We are using the following set of flags:
* GCC:
   ````
      -Wall -Wunused-parameter -Wshadow -pedantic
      -Werror -Wsign-compare -Wtype-limits
      -Wignored-qualifiers -Wempty-body
      -Wclobbered -Wuninitialized
   ````
* Intel:
    ````
      -Wall -Wunused-parameter -Wshadow -pedantic
      -Werror -Wsign-compare -Wtype-limits
      -Wuninitialized
    ````
* Clang:
    ````
      -Wall -Wunused-parameter -Wshadow -pedantic
      -Werror -Wsign-compare -Wtype-limits
      -Wuninitialized
    ````

* NVCC:
  ````
    -Wall -Wunused-parameter -Wshadow -pedantic
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
````bash
cmake $srcdir \
  -DCMAKE_CXX_COMPILER=$path_to_compiler \
  -DCMAKE_INSTALL_PREFIX=$path_to_install \
  -DKokkos_ENABLE_OPENMP=On \
  -DKokkos_ARCH_HSW=On \
  -DKokkos_HWLOC_DIR=$path_to_hwloc
````
then simply type `make install`. The Kokkos CMake package will then be installed in `$path_to_install` to be used by downstream packages.

To validate the Kokkos build, configure with
````
 -DKokkos_ENABLE_TESTS=On
````
and run `make test` after completing the build.

For your CMake project using Kokkos, code such as the following:

````cmake
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
To get started, download the Spack [repo](https://github.com/spack/spack).
````
A basic installation would be done as:
````bash
> spack install kokkos
````
Spack allows options and and compilers to be tuned in the install command.
````bash
> spack install kokkos@3.0 %gcc@7.3.0 +openmp
````
This example illustrates the three most common parameters to Spack:
* Variants: specified with, e.g. `+openmp`, this activates (or deactivates with, e.g. `~openmp`) certain options.
* Version:  immediately following `kokkos` the `@version` can specify a particular Kokkos to build
* Compiler: a default compiler will be chosen if not specified, but an exact compiler version can be given with the `%`option.

For a complete list of Kokkos options, run:
````bash
> spack info kokkos
````
Spack currently installs packages to a location determined by a unique hash. This hash name is not really "human readable".
Generally, Spack usage should never really require you to reference the computer-generated unique install folder.
More details are given in the [build instructions](BUILD.md). If you must know, you can locate Spack Kokkos installations with:
````bash
> spack find -p kokkos ...
````
where `...` is the unique spec identifying the particular Kokkos configuration and version.
Some more details can found in the Kokkos spack [documentation](Spack.md) or the Spack [website](https://spack.readthedocs.io/en/latest).

## Raw Makefile

Raw Makefiles are only supported via inline builds. See below.

## Inline Builds vs. Installed Package
For individual projects, it may be preferable to build Kokkos inline rather than link to an installed package.
The main reason is that you may otherwise need many different
configurations of Kokkos installed depending on the required compile time
features an application needs. For example there is only one default
execution space, which means you need different installations to have OpenMP
or C++ threads as the default space. Also for the CUDA backend there are certain
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

````BibTex
@ARTICLE{9485033,
  author={Trott, Christian R. and Lebrun-Grandié, Damien and Arndt, Daniel and Ciesko, Jan and Dang, Vinh and Ellingwood, Nathan and Gayatri, Rahulkumar and Harvey, Evan and Hollman, Daisy S. and Ibanez, Dan and Liber, Nevin and Madsen, Jonathan and Miles, Jeff and Poliakoff, David and Powell, Amy and Rajamanickam, Sivasankaran and Simberg, Mikael and Sunderland, Dan and Turcksin, Bruno and Wilke, Jeremiah},
  journal={IEEE Transactions on Parallel and Distributed Systems},
  title={Kokkos 3: Programming Model Extensions for the Exascale Era},
  year={2022},
  volume={33},
  number={4},
  pages={805-817},
  doi={10.1109/TPDS.2021.3097283}}
````

If you use more than one Kokkos EcoSystem package, please also cite:

````BibTex
@ARTICLE{9502936,
  author={Trott, Christian and Berger-Vergiat, Luc and Poliakoff, David and Rajamanickam, Sivasankaran and Lebrun-Grandie, Damien and Madsen, Jonathan and Al Awar, Nader and Gligoric, Milos and Shipman, Galen and Womeldorff, Geoff},
  journal={Computing in Science   Engineering},
  title={The Kokkos EcoSystem: Comprehensive Performance Portability for High Performance Computing},
  year={2021},
  volume={23},
  number={5},
  pages={10-18},
  doi={10.1109/MCSE.2021.3098509}}
````


And if you feel generous: feel free to cite the original Kokkos paper which describes most of the basic Kokkos concepts:

````BibTeX
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

