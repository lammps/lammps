Kokkos Core implements a programming model in C++ for writing performance portable
applications targeting all major HPC platforms. For that purpose it provides
abstractions for both parallel execution of code and data management.
Kokkos is designed to target complex node architectures with N-level memory
hierarchies and multiple types of execution resources. It currently can use
OpenMP, Pthreads and CUDA as backend programming models.

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
  * PGI 18.7
  * NVCC 7.5 for CUDA (with gcc 4.8.4)
  * NVCC 8.0.44 for CUDA (with gcc 5.3.0)
  * NVCC 9.1 for CUDA (with gcc 6.1.0)
  * NVCC 9.2 for CUDA (with gcc 7.2.0)
  * NVCC 10.0 for CUDA (with gcc 7.4.0)

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
   - Cygwin 2.1.0 64bit with gcc 4.9.3
   - GCC 8.1.0 (not warning free)

### Known non-working combinations:
  * Power8:
   - Pthreads backend
  * ARM
   - Pthreads backend


Primary tested compiler are passing in release mode
with warnings as errors. They also are tested with a comprehensive set of 
backend combinations (i.e. OpenMP, Pthreads, Serial, OpenMP+Serial, ...).
We are using the following set of flags:
GCC:   -Wall -Wshadow -pedantic -Werror -Wsign-compare -Wtype-limits
       -Wignored-qualifiers -Wempty-body -Wclobbered -Wuninitialized
Intel: -Wall -Wshadow -pedantic -Werror -Wsign-compare -Wtype-limits -Wuninitialized
Clang: -Wall -Wshadow -pedantic -Werror -Wsign-compare -Wtype-limits -Wuninitialized
NVCC:  -Wall -Wshadow -pedantic -Werror -Wsign-compare -Wtype-limits -Wuninitialized

Other compilers are tested occasionally, in particular when pushing from develop to 
master branch, without -Werror and only for a select set of backends.

# Running Unit Tests

To run the unit tests create a build directory and run the following commands

KOKKOS_PATH/generate_makefile.bash
make build-test
make test

Run KOKKOS_PATH/generate_makefile.bash --help for more detailed options such as
changing the device type for which to build.

# Installing the library

To install Kokkos as a library create a build directory and run the following

KOKKOS_PATH/generate_makefile.bash --prefix=INSTALL_PATH
make kokkoslib
make install

KOKKOS_PATH/generate_makefile.bash --help for more detailed options such as
changing the device type for which to build.

Note that in many cases it is preferable to build Kokkos inline with an 
application. The main reason is that you may otherwise need many different
configurations of Kokkos installed depending on the required compile time
features an application needs. For example there is only one default 
execution space, which means you need different installations to have OpenMP
or Pthreads as the default space. Also for the CUDA backend there are certain
choices, such as allowing relocatable device code, which must be made at 
installation time. Building Kokkos inline uses largely the same process
as compiling an application against an installed Kokkos library. See for 
example benchmarks/bytes_and_flops/Makefile which can be used with an installed
library and for an inline build.  

### CMake

Kokkos supports being build as part of a CMake applications. An example can 
be found in example/cmake_build. 

# Kokkos and CUDA UVM

Kokkos does support UVM as a specific memory space called CudaUVMSpace. 
Allocations made with that space are accessible from host and device. 
You can tell Kokkos to use that as the default space for Cuda allocations.
In either case UVM comes with a number of restrictions:
(i) You can't access allocations on the host while a kernel is potentially 
running. This will lead to segfaults. To avoid that you either need to 
call Kokkos::Cuda::fence() (or just Kokkos::fence()), after kernels, or
you can set the environment variable CUDA_LAUNCH_BLOCKING=1.
Furthermore in multi socket multi GPU machines without NVLINK, UVM defaults 
to using zero copy allocations for technical reasons related to using multiple
GPUs from the same process. If an executable doesn't do that (e.g. each
MPI rank of an application uses a single GPU [can be the same GPU for 
multiple MPI ranks]) you can set CUDA_MANAGED_FORCE_DEVICE_ALLOC=1.
This will enforce proper UVM allocations, but can lead to errors if 
more than a single GPU is used by a single process.


# Citing Kokkos

If you publish work which mentions Kokkos, please cite the following paper:

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
