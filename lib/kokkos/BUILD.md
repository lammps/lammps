![Kokkos](https://avatars2.githubusercontent.com/u/10199860?s=200&v=4)

# Installing and Using Kokkos

## Kokkos Philosophy
Kokkos provides a modern CMake style build system.
As C++ continues to develop for C++20 and beyond, CMake is likely to provide the most robust support
for C++.  Applications heavily leveraging Kokkos are strongly encouraged to use a CMake build system.

You can either use Kokkos as an installed package (encouraged) or use Kokkos in-tree in your project.
Modern CMake is exceedingly simple at a high-level (with the devil in the details).
Once Kokkos is installed In your `CMakeLists.txt` simply use:
````cmake
find_package(Kokkos REQUIRED)
````
Then for every executable or library in your project:
````cmake
target_link_libraries(myTarget Kokkos::kokkos)
````
That's it! There is no checking Kokkos preprocessor, compiler, or linker flags.
Kokkos propagates all the necessary flags to your project.
This means not only is linking to Kokkos easy, but Kokkos itself can actually configure compiler and linker flags for *your*
project.
When configuring your project just set:
````bash
> cmake ${srcdir} \
  -DKokkos_ROOT=${kokkos_install_prefix} \
  -DCMAKE_CXX_COMPILER=${compiler_used_to_build_kokkos}
````
Note: You may need the following if using some versions of CMake (e.g. 3.12):
````cmake
cmake_policy(SET CMP0074 NEW)
````
If building in-tree, there is no `find_package`. You can use `add_subdirectory(kokkos)` with the Kokkos source and again just link with `target_link_libraries(Kokkos::kokkos)`.
The examples in `examples/cmake_build_installed` and `examples/cmake_build_in_tree` can help get you started.


## Configuring CMake
A very basic installation of Kokkos is done with:
````bash
> cmake ${srcdir} \
 -DCMAKE_CXX_COMPILER=g++ \
 -DCMAKE_INSTALL_PREFIX=${kokkos_install_folder}
````
which builds and installed a default Kokkos when you run `make install`.
There are numerous device backends, options, and architecture-specific optimizations that can be configured, e.g.
````bash
> cmake ${srcdir} \
 -DCMAKE_CXX_COMPILER=g++ \
 -DCMAKE_INSTALL_PREFIX=${kokkos_install_folder} \
 -DKokkos_ENABLE_OPENMP=ON
````
which activates the OpenMP backend. All of the options controlling device backends, options, architectures, and third-party libraries (TPLs) are given below.

## Known Issues<a name="KnownIssues"></a>

### Cray

* The Cray compiler wrappers do static linking by default. This seems to break the Kokkos build. You will likely need to set the environment variable `CRAYPE_LINK_TYPE=dynamic` in order to link correctly. Kokkos warns during configure if this is missing.
* The Cray compiler identifies to CMake as Clang, but it sometimes has its own flags that differ from Clang. We try to include all exceptions, but flag errors may occur in which a Clang-specific flag is passed that the Cray compiler does not recognize.

### Fortran

* In a mixed C++/Fortran code, CMake will use the C++ linker by default. If you override this behavior and use Fortran as the link language, the link may break because Kokkos adds linker flags expecting the linker to be C++. Prior to CMake 3.18, Kokkos has no way of detecting in downstream projects that the linker was changed to Fortran.  From CMake 3.18, Kokkos can use generator expressions to avoid adding flags when the linker is not C++. Note: Kokkos will not add any linker flags in this Fortran case. The user will be entirely on their own to add the appropriate linker flags.

## Spack
An alternative to manually building with the CMake is to use the Spack package manager.
Make sure you have downloaded [Spack](https://github.com/spack/spack).
The easiest way to configure the Spack environment is:
````bash
> source spack/share/spack/setup-env.sh
````
with other scripts available for other shells.
You can display information about how to install packages with:
````bash
> spack info kokkos
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
More details can be found in the [Spack README](Spack.md)

#### Spack Development
Spack currently installs packages to a location determined by a unique hash. This hash name is not really "human readable".
Generally, Spack usage should never really require you to reference the computer-generated unique install folder.
If you must know, you can locate Spack Kokkos installations with:
````bash
> spack find -p kokkos ...
````
where `...` is the unique spec identifying the particular Kokkos configuration and version.

A better way to use Spack for doing Kokkos development is the dev-build feature of Spack.
For dev-build details, consult the kokkos-spack repository [README](https://github.com/kokkos/kokkos-spack/blob/master/README.md).

# Kokkos Keyword Listing

## Device Backends
Device backends can be enabled by specifying `-DKokkos_ENABLE_X`.

* Kokkos_ENABLE_CUDA
    * Whether to build CUDA backend
    * BOOL Default: OFF
* Kokkos_ENABLE_HPX
    * Whether to build HPX backend (experimental)
    * BOOL Default: OFF
* Kokkos_ENABLE_OPENMP
    * Whether to build OpenMP backend
    * BOOL Default: OFF
* Kokkos_ENABLE_PTHREAD
    * Whether to build Pthread backend
    * BOOL Default: OFF
* Kokkos_ENABLE_SERIAL
    * Whether to build serial backend
    * BOOL Default: ON
* Kokkos_ENABLE_HIP (Experimental)
    * Whether to build HIP backend
    * BOOL Default: OFF
* Kokkos_ENABLE_OPENMPTARGET (Experimental)
    * Whether to build the OpenMP target backend
    * BOOL Default: OFF

## Enable Options
Options can be enabled by specifying `-DKokkos_ENABLE_X`.

* Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION
    * Whether to aggressively vectorize loops
    * BOOL Default: OFF
* Kokkos_ENABLE_COMPILER_WARNINGS
    * Whether to print all compiler warnings
    * BOOL Default: OFF
* Kokkos_ENABLE_CUDA_CONSTEXPR
    * Whether to activate experimental relaxed constexpr functions
    * BOOL Default: OFF
* Kokkos_ENABLE_CUDA_LAMBDA
    * Whether to activate experimental lambda features
    * BOOL Default: OFF
* Kokkos_ENABLE_CUDA_LDG_INTRINSIC
    * Whether to use CUDA LDG intrinsics
    * BOOL Default: OFF
* Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE
    * Whether to enable relocatable device code (RDC) for CUDA
    * BOOL Default: OFF
* Kokkos_ENABLE_CUDA_UVM
    * Whether to use unified memory (UM) by default for CUDA
    * BOOL Default: OFF
* Kokkos_ENABLE_DEBUG
    * Whether to activate extra debug features - may increase compile times
    * BOOL Default: OFF
* Kokkos_ENABLE_DEBUG_BOUNDS_CHECK
    * Whether to use bounds checking - will increase runtime
    * BOOL Default: OFF
* Kokkos_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
    * Debug check on dual views
    * BOOL Default: OFF
* Kokkos_ENABLE_EXAMPLES
    * Whether to enable building examples
    * BOOL Default: OFF
* Kokkos_ENABLE_HPX_ASYNC_DISPATCH
    * Whether HPX supports asynchronous dispatch
    * BOOL Default: OFF
* Kokkos_ENABLE_LARGE_MEM_TESTS
    * Whether to perform extra large memory tests
    * BOOL_Default: OFF
* Kokkos_ENABLE_PROFILING_LOAD_PRINT
    * Whether to print information about which profiling tools gotloaded
    * BOOL Default: OFF
* Kokkos_ENABLE_TESTS
    * Whether to build serial  backend
    * BOOL Default: OFF

## Other Options
* Kokkos_CXX_STANDARD
    * The C++ standard for Kokkos to use: c++14, c++17, or c++20. This should be given in CMake style as 14, 17, or 20.
    * STRING Default: 14

## Third-party Libraries (TPLs)
The following options control enabling TPLs:
* Kokkos_ENABLE_HPX
    * Whether to enable the HPX library
    * BOOL Default: OFF
* Kokkos_ENABLE_HWLOC
    * Whether to enable the HWLOC library
    * BOOL Default: Off
* Kokkos_ENABLE_LIBNUMA
    * Whether to enable the LIBNUMA library
    * BOOL Default: Off
* Kokkos_ENABLE_MEMKIND
    * Whether to enable the MEMKIND library
    * BOOL Default: Off
* Kokkos_ENABLE_LIBDL
    * Whether to enable the LIBDL library
    * BOOL Default: On
* Kokkos_ENABLE_LIBRT
    * Whether to enable the LIBRT library
    * BOOL Default: Off

The following options control finding and configuring non-CMake TPLs:
* Kokkos_CUDA_DIR or CUDA_ROOT
    * Location of CUDA install prefix for libraries
    * PATH Default:
* Kokkos_HWLOC_DIR or HWLOC_ROOT
    * Location of HWLOC install prefix
    * PATH Default:
* Kokkos_LIBNUMA_DIR or LIBNUMA_ROOT
    * Location of LIBNUMA install prefix
    * PATH Default:
* Kokkos_MEMKIND_DIR or MEMKIND_ROOT
    * Location of MEMKIND install prefix
    * PATH Default:
* Kokkos_LIBDL_DIR or LIBDL_ROOT
    * Location of LIBDL install prefix
    * PATH Default:
* Kokkos_LIBRT_DIR or LIBRT_ROOT
    * Location of LIBRT install prefix
    * PATH Default:

The following options control `find_package` paths for CMake-based TPLs:
* HPX_DIR or HPX_ROOT
    * Location of HPX prefix (ROOT) or CMake config file (DIR)
    * PATH Default:

## Architecture Keywords
Architecture-specific optimizations can be enabled by specifying `-DKokkos_ARCH_X`.

* Kokkos_ARCH_AMDAVX
    * Whether to optimize for the AMDAVX architecture
    * BOOL Default: OFF
* Kokkos_ARCH_ARMV80
    * Whether to optimize for the ARMV80 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_ARMV81
    * Whether to optimize for the ARMV81 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_ARMV8_THUNDERX
    * Whether to optimize for the ARMV8_THUNDERX architecture
    * BOOL Default: OFF
* Kokkos_ARCH_ARMV8_TX2
    * Whether to optimize for the ARMV8_TX2 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_BDW
    * Whether to optimize for the BDW architecture
    * BOOL Default: OFF
* Kokkos_ARCH_BGQ
    * Whether to optimize for the BGQ architecture
    * BOOL Default: OFF
* Kokkos_ARCH_ZEN
    * Whether to optimize for the Zen architecture
    * BOOL Default: OFF
* Kokkos_ARCH_ZEN2
    * Whether to optimize for the Zen2 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_HSW
    * Whether to optimize for the HSW architecture
    * BOOL Default: OFF
* Kokkos_ARCH_KEPLER30
    * Whether to optimize for the KEPLER30 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_KEPLER32
    * Whether to optimize for the KEPLER32 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_KEPLER35
    * Whether to optimize for the KEPLER35 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_KEPLER37
    * Whether to optimize for the KEPLER37 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_KNC
    * Whether to optimize for the KNC architecture
    * BOOL Default: OFF
* Kokkos_ARCH_KNL
    * Whether to optimize for the KNL architecture
    * BOOL Default: OFF
* Kokkos_ARCH_MAXWELL50
    * Whether to optimize for the MAXWELL50 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_MAXWELL52
    * Whether to optimize for the MAXWELL52 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_MAXWELL53
    * Whether to optimize for the MAXWELL53 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_PASCAL60
    * Whether to optimize for the PASCAL60 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_PASCAL61
    * Whether to optimize for the PASCAL61 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_POWER7
    * Whether to optimize for the POWER7 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_POWER8
    * Whether to optimize for the POWER8 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_POWER9
    * Whether to optimize for the POWER9 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_SKX
    * Whether to optimize for the SKX architecture
    * BOOL Default: OFF
* Kokkos_ARCH_SNB
    * Whether to optimize for the SNB architecture
    * BOOL Default: OFF
* Kokkos_ARCH_TURING75
    * Whether to optimize for the TURING75 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_VOLTA70
    * Whether to optimize for the VOLTA70 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_VOLTA72
    * Whether to optimize for the VOLTA72 architecture
    * BOOL Default: OFF
* Kokkos_ARCH_WSM
    * Whether to optimize for the WSM architecture
    * BOOL Default: OFF

##### [LICENSE](https://github.com/kokkos/kokkos/blob/devel/LICENSE)

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Under the terms of Contract DE-NA0003525 with NTESS,
the U.S. Government retains certain rights in this software.
