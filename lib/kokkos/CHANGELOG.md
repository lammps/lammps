# CHANGELOG

## [4.3.00](https://github.com/kokkos/kokkos/tree/4.3.00) (2024-03-19)
[Full Changelog](https://github.com/kokkos/kokkos/compare/4.2.01...4.3.00)

### Features:
* Add `Experimental::sort_by_key(exec, keys, values)` algorithm [\#6801](https://github.com/kokkos/kokkos/pull/6801)

### Backend and Architecture Enhancements:

#### CUDA:
* Experimental multi-GPU support (from the same process) [\#6782](https://github.com/kokkos/kokkos/pull/6782)
* Link against CUDA libraries even with KOKKOS_ENABLE_COMPILE_AS_CMAKE_LANGUAGE [\#6701](https://github.com/kokkos/kokkos/pull/6701)
* Don't use the compiler launcher script if the CMake compile language is CUDA. [\#6704](https://github.com/kokkos/kokkos/pull/6704)
* nvcc(wrapper): adding "long" and "short" versions for all flags [\#6615](https://github.com/kokkos/kokkos/pull/6615)

#### HIP:
 * Fix compilation when using amdclang (with ROCm >= 5.7) and RDC [\#6857](https://github.com/kokkos/kokkos/pull/6857)
 * Use rocthrust for sorting, when available [\#6793](https://github.com/kokkos/kokkos/pull/6793)

#### SYCL:
* We only support OneAPI SYCL implementation: add check during initialization
  * Error out on initialization if the backend is different from `ext_oneapi_*` [\#6784](https://github.com/kokkos/kokkos/pull/6784)
  * Filter GPU devices for `ext_onapi_*` GPU devices [\#6758](https://github.com/kokkos/kokkos/pull/6784)
* Performance Improvements
  * Avoid unnecessary zero-memset of the scratch flags in SYCL [\#6739](https://github.com/kokkos/kokkos/pull/6739)
  * Use host-pinned memory to copy reduction/scan result [\#6500](https://github.com/kokkos/kokkos/pull/6500)
* Address deprecations after oneAPI 2023.2.0 [\#6577](https://github.com/kokkos/kokkos/pull/6739)
* Make sure to call find_dependency for oneDPL if necessary [\#6870](https://github.com/kokkos/kokkos/pull/6870)

#### OpenMPTarget:
* Use LLVM extensions for dynamic shared memory [\#6380](https://github.com/kokkos/kokkos/pull/6380)
* Guard scratch memory usage in ParallelReduce [\#6585 ](https://github.com/kokkos/kokkos/pull/6585)
* Update linker flags for Intel GPUs update [\#6735](https://github.com/kokkos/kokkos/pull/6735)
* Improve handling of printf on Intel GPUs [\#6652](https://github.com/kokkos/kokkos/pull/6652)

#### OpenACC:
* Add atomics support [\#6446](https://github.com/kokkos/kokkos/pull/6446)
* Make the OpenACC backend asynchronous [\#6772](https://github.com/kokkos/kokkos/pull/6772)

#### Threads:
* Add missing broadcast to TeamThreadRange parallel_scan [\#6601](https://github.com/kokkos/kokkos/pull/6446)

#### OpenMP:
* Improve performance of view initializations and filling with zeros [\#6573](https://github.com/kokkos/kokkos/pull/6573)

### General Enhancements

* Improve performance of random number generation when using a normal distribution on GPUs [\#6556](https://github.com/kokkos/kokkos/pull/6556)
* Allocate temporary view with the user-provided execution space instance and do not initialize in `unique` algorithm [\#6598](https://github.com/kokkos/kokkos/pull/6598)
* Add deduction guide for `Kokkos::Array` [\#6373](https://github.com/kokkos/kokkos/pull/6373)
* Provide new public headers `<Kokkos_Clamp.hpp>` and `<Kokkos_MinMax.hpp>` [\#6687](https://github.com/kokkos/kokkos/pull/6687)
* Fix/improvement to `remove_if` parallel algorithm: use the provided execution space instance for temporary allocations and drop unnecessaryinitialization + avoid evaluating twice the predicate during final pass [\#6747](https://github.com/kokkos/kokkos/pull/6747)
* Add runtime function to query the number of devices and make device ID consistent with `KOKKOS_VISIBLE_DEVICES` [\#6713](https://github.com/kokkos/kokkos/pull/6713)
* simd: support `vector_aligned_tag` [\#6243](https://github.com/kokkos/kokkos/pull/6243)
* Avoid unnecessary allocation when default constructing Bitset [\#6524](https://github.com/kokkos/kokkos/pull/6524)
* Fix constness for views in std algorithms [\#6813](https://github.com/kokkos/kokkos/pull/6813)
* Improve error message on unsafe implicit conversion in MDRangePolicy [\#6855](https://github.com/kokkos/kokkos/pull/6855)
* CTAD (deduction guides) for RangePolicy [\#6850](https://github.com/kokkos/kokkos/pull/6850)
* CTAD (deduction guides) for MDRangePolicy [\#5516](https://github.com/kokkos/kokkos/pull/5516)

### Build System Changes
* Require `Kokkos_ENABLE_ATOMICS_BYPASS` option to bypass atomic operation for Serial backend only builds [\#6692](https://github.com/kokkos/kokkos/pull/6692)
* Add support for RISCV and the Milk-V's Pioneer [\#6773](https://github.com/kokkos/kokkos/pull/6773)
* Add C++26 standard to CMake setup [\#6733](https://github.com/kokkos/kokkos/pull/6733)
* Fix Makefile when using gnu_generate_makefile.sh and make >= 4.3 [\#6606](https://github.com/kokkos/kokkos/pull/6606)
* Cuda: Fix configuring with CMake >= 3.28.4 - temporary fallback to internal CudaToolkit.cmake [\#6898](https://github.com/kokkos/kokkos/pull/6898)

### Incompatibilities (i.e. breaking changes)
* Remove all `DEPRECATED_CODE_3` option and all code that was guarded by it  [\#6523](https://github.com/kokkos/kokkos/pull/6523)
* Drop guards to accommodate external code defining `KOKKOS_ASSERT` [\#6665](https://github.com/kokkos/kokkos/pull/6665)
* `Profiling::ProfilingSection(std::string)` constructor marked explicit and nodiscard [\#6690](https://github.com/kokkos/kokkos/pull/6690)
* Add bound check preconditions for `RangePolicy` and `MDRangePolicy` [\#6617](https://github.com/kokkos/kokkos/pull/6617) [\#6726](https://github.com/kokkos/kokkos/pull/6726)
* Add checks for unsafe implicit conversions in RangePolicy [\#6754](https://github.com/kokkos/kokkos/pull/6754)
* Remove Kokkos::[b]half_t volatile overloads [\#6579](https://github.com/kokkos/kokkos/pull/6579)
* Remove KOKKOS_IMPL_DO_NOT_USE_PRINTF [\#6593](https://github.com/kokkos/kokkos/pull/6593)
* Check matching static extents in View constructor [\#5190 ](https://github.com/kokkos/kokkos/pull/5190)
* Tools(profiling): fix typo Kokkos_Tools_Optim[i]zationGoal [\#6642](https://github.com/kokkos/kokkos/pull/6642)
* Remove variadic range policy constructor (disallow passing multiple trailing chunk size arguments) [\#6845](https://github.com/kokkos/kokkos/pull/6845)
* Improve message on view out of bounds access and always abort [\#6861](https://github.com/kokkos/kokkos/pull/6861)
* Drop `KOKKOS_ENABLE_INTEL_MM_ALLOC` macro [\#6797](https://github.com/kokkos/kokkos/pull/6797)
* Remove `Kokkos::Experimental::LogicalMemorySpace` (without going through deprecation) [\#6557](https://github.com/kokkos/kokkos/pull/6557)
* Remove `Experimental::HBWSpace` and support for linking against memkind [\#6791](https://github.com/kokkos/kokkos/pull/6791)
* Drop librt TPL and associated `KOKKOS_ENABLE_LIBRT` macro [\#6798](https://github.com/kokkos/kokkos/pull/6798)
* Drop support for old CPU architectures (`ARCH_BGQ`, `ARCH_POWER7`, `ARCH_WSM` and associated `ARCH_SSE4` macro) [\#6806](https://github.com/kokkos/kokkos/pull/6806)
* Drop support for deprecated command-line arguments and environment variables [\#6744](https://github.com/kokkos/kokkos/pull/6744)

### Deprecations
* Provide kokkos_swap as part of Core and deprecate Experimental::swap in Algorithms [\#6697](https://github.com/kokkos/kokkos/pull/6697)
* Deprecate {Cuda,HIP}::detect_device_count() and Cuda::[detect_]device_arch() [\#6710](https://github.com/kokkos/kokkos/pull/6710)
* Deprecate `ExecutionSpace::in_parallel()` [\#6582](https://github.com/kokkos/kokkos/pull/6582)

### Bug Fixes
* Fix team-level MDRange reductions: [\#6511](https://github.com/kokkos/kokkos/pull/6511)
* Fix CUDA and SYCL small value type (16-bit) team reductions [\#5334](https://github.com/kokkos/kokkos/pull/5334)
* Enable `{transform_}exclusive_scan` in place [\#6667](https://github.com/kokkos/kokkos/pull/6667)
* `fill_random` overload that do not take an execution space instance argument should fence [\#6658](https://github.com/kokkos/kokkos/pull/6658)
* HIP,Cuda,OpenMPTarget: Fixup use provided execution space when copying host inaccessible reduction result [\#6777](https://github.com/kokkos/kokkos/pull/6777)
* Fix typo in `cuda_func_set_attribute[s]_wrapper` preventing proper setting of desired occupancy [\#6786](https://github.com/kokkos/kokkos/pull/6786)
* Avoid undefined behavior due to conversion between signed and unsigned integers in shift_{right, left}_team_impl [\#6821](https://github.com/kokkos/kokkos/pull/6821)
* Fix a bug in Makefile.kokkos when using AMD GPU architectures as `AMD_GFXYYY` [\#6892](https://github.com/kokkos/kokkos/pull/6892)

## [4.2.01](https://github.com/kokkos/kokkos/tree/4.2.01) (2023-12-07)
[Full Changelog](https://github.com/kokkos/kokkos/compare/4.2.00...4.2.01)

### Backend and Architecture Enhancements:

#### CUDA:
- Add warp sync for `parallel_reduce` to avoid race condition [\#6630](https://github.com/kokkos/kokkos/pull/6630), [\#6746](https://github.com/kokkos/kokkos/pull/6746)

#### HIP:
- Fix Graph "multiple definition of" linking error (missing `inline` specifier) [\#6624](https://github.com/kokkos/kokkos/pull/6624)
- Add support for gfx940 (AMD Instinct MI300 GPU) [\#6671](https://github.com/kokkos/kokkos/pull/6671)

### Build System
- CMake: Don't let Kokkos set `CMAKE_CXX_FLAGS` for Trilinos builds [\#6742](https://github.com/kokkos/kokkos/pull/6742)

### Bug Fixes
- Remove deprecation warning for `AllocationMechanism` for GCC <11.0 [\#6653](https://github.com/kokkos/kokkos/pull/6653)
- Fix bug early tools finalize with non-default host execution instances [\#6635](https://github.com/kokkos/kokkos/pull/6635)
- Fix various issues for MSVC CUDA builds [\#6659](https://github.com/kokkos/kokkos/pull/6659)
- Fix "extra `;`" warning with `-pedantic` flag in `<Kokkos_SIMD_Scalar.hpp>` [\#6510](https://github.com/kokkos/kokkos/pull/6510)

## [4.2.00](https://github.com/kokkos/kokkos/tree/4.2.00) (2023-11-06)
[Full Changelog](https://github.com/kokkos/kokkos/compare/4.1.00...4.2.00)

### Features:
- SIMD: significant improvements to SIMD support and alignment with C++26 SIMD
  - add `Kokkos::abs` overload for SIMD types [\#6069](https://github.com/kokkos/kokkos/pull/6069)
  - add generator constructors [\#6347](https://github.com/kokkos/kokkos/pull/6347)
  - convert binary operators to hidden friends [\#6320](https://github.com/kokkos/kokkos/pull/6320)
  - add shift operators [\#6109](https://github.com/kokkos/kokkos/pull/6109)
  - add `float` support [\#6177](https://github.com/kokkos/kokkos/pull/6177)
  - add remaining `gather_from` and `scatter_to` overloads [\#6220](https://github.com/kokkos/kokkos/pull/6220)
  - define simd math function overloads in the Kokkos namespace [\#6465](https://github.com/kokkos/kokkos/pull/6465), [\#6487](https://github.com/kokkos/kokkos/pull/6487)
  - `Kokkos_ENABLE_NATIVE=ON` autodetects SIMD types supported [\#6188](https://github.com/kokkos/kokkos/pull/6188)
  - fix AVX2 SIMD support for ZEN2 AMD CPU [\#6238](https://github.com/kokkos/kokkos/pull/6238)
- `Kokkos::printf` [\#6083](https://github.com/kokkos/kokkos/pull/6083)
- `Kokkos::sort`: support custom comparator [\#6253](https://github.com/kokkos/kokkos/pull/6253)
- `half_t` and `bhalf_t` numeric traits [\#5778](https://github.com/kokkos/kokkos/pull/5778)
- `half_t` and `bhalf_t` mixed comparisons [\#6407](https://github.com/kokkos/kokkos/pull/6407)
- `half_t` and `bhalf_t` mathematical functions [\#6124](https://github.com/kokkos/kokkos/pull/6124)
- `TeamThreadRange` `parallel_scan` with return value [\#6090](https://github.com/kokkos/kokkos/pull/6090), [\#6301](https://github.com/kokkos/kokkos/pull/6301), [\#6302](https://github.com/kokkos/kokkos/pull/6302), [\#6303](https://github.com/kokkos/kokkos/pull/6303), [\#6307](https://github.com/kokkos/kokkos/pull/6307)
- `ThreadVectorRange` `parallel_scan` with return value [\#6235](https://github.com/kokkos/kokkos/pull/6235), [\#6242](https://github.com/kokkos/kokkos/pull/6242), [\#6308](https://github.com/kokkos/kokkos/pull/6308), [\#6305](https://github.com/kokkos/kokkos/pull/6305), [\#6292](https://github.com/kokkos/kokkos/pull/6292)
- Add team-level std algorithms [\#6200](https://github.com/kokkos/kokkos/pull/6200), [\#6205](https://github.com/kokkos/kokkos/pull/6205), [\#6207](https://github.com/kokkos/kokkos/pull/6207), [\#6208](https://github.com/kokkos/kokkos/pull/6208), [\#6209](https://github.com/kokkos/kokkos/pull/6209), [\#6210](https://github.com/kokkos/kokkos/pull/6210), [\#6211](https://github.com/kokkos/kokkos/pull/6211), [\#6212](https://github.com/kokkos/kokkos/pull/6212), [\#6213](https://github.com/kokkos/kokkos/pull/6213), [\#6256](https://github.com/kokkos/kokkos/pull/6256), [\#6258](https://github.com/kokkos/kokkos/pull/6258), [\#6350](https://github.com/kokkos/kokkos/pull/6350), [\#6351](https://github.com/kokkos/kokkos/pull/6351)
- Serial: Allow for distinct execution space instances [\#6441](https://github.com/kokkos/kokkos/pull/6441)

### Backend and Architecture Enhancements:

#### CUDA:
- Fixed potential data race in Cuda `parallel_reduce` [\#6236](https://github.com/kokkos/kokkos/pull/6236)
- Use `cudaMallocAsync` by default [\#6402](https://github.com/kokkos/kokkos/pull/6402)
- Bugfix for using Kokkos from a thread of execution [\#6299](https://github.com/kokkos/kokkos/pull/6299)

#### HIP:
- New naming convention for AMD GPU: VEGA906, VEGA908, VEGA90A, NAVI1030 to AMD_GFX906, AMD_GFX908, AMD_GFX90A, AMD_GFX1030 [\#6266](https://github.com/kokkos/kokkos/pull/6266)
- Add initial support for gfx942: [\#6358](https://github.com/kokkos/kokkos/pull/6358)
- Improve reduction performance [\#6229](https://github.com/kokkos/kokkos/pull/6229)
- Deprecate `HIP(hipStream_t,bool)` constructor [\#6401](https://github.com/kokkos/kokkos/pull/6401)
- Add support for Graph [\#6370](https://github.com/kokkos/kokkos/pull/6370)
- Improve reduction performance when using Teams [\#6284](https://github.com/kokkos/kokkos/pull/6284)
- Fix concurrency calculation [\#6479](https://github.com/kokkos/kokkos/pull/6479)
- Fix potential data race in HIP `parallel_reduce` [\#6429](https://github.com/kokkos/kokkos/pull/6429)

#### SYCL:
- Enforce external `sycl::queues` to be in-order [\#6246](https://github.com/kokkos/kokkos/pull/6246)
- Improve reduction performance: [\#6272](https://github.com/kokkos/kokkos/pull/6272) [\#6271](https://github.com/kokkos/kokkos/pull/6271) [\#6270](https://github.com/kokkos/kokkos/pull/6270) [\#6264](https://github.com/kokkos/kokkos/pull/6264)
- Allow using the SYCL execution space on AMD GPUs [\#6321](https://github.com/kokkos/kokkos/pull/6321)
- Allow sorting via native oneDPL to support Views with stride=1 [\#6322](https://github.com/kokkos/kokkos/pull/6322)
- Make in-order queues the default via macro [\#6189](https://github.com/kokkos/kokkos/pull/6189)

#### OpenACC:
- Support Clacc compiler [\#6250](https://github.com/kokkos/kokkos/pull/6250)

### General Enhancements
- Add missing `is_*_view` traits and `is_*_view_v` helper variable templates for `DynRankView`, `DynamicView`, `OffsetView`, `ScatterView` containers [\#6195](https://github.com/kokkos/kokkos/pull/6195)
- Make `nvcc_wrapper` and `compiler_launcher` scripts more portable by switching to a `#!/usr/bin/env` shebang [\#6357](https://github.com/kokkos/kokkos/pull/6357)
- Add an improved `Kokkos::malloc` / `Kokkos::free` performance test [\#6377](https://github.com/kokkos/kokkos/pull/6377)
- Ensure `Views` with `size==0` can be used with `deep_copy` [\#6273](https://github.com/kokkos/kokkos/pull/6273)
- `Kokkos::abort` is moved to header `Kokkos_Abort.hpp` [\#6445](https://github.com/kokkos/kokkos/pull/6445)
- `KOKKOS_ASSERT`, `KOKKOS_EXPECTS`, `KOKKOS_ENSURES` are moved to header `Kokkos_Assert.hpp` [\#6445](https://github.com/kokkos/kokkos/pull/6445)
- Add a permuted-index mode to the gups benchmark [\#6378](https://github.com/kokkos/kokkos/pull/6378)
- Check for overflow during backend initialization [\#6159](https://github.com/kokkos/kokkos/pull/6159)
- Make constraints on `Kokkos::sort` more visible [\#6234](https://github.com/kokkos/kokkos/pull/6234) and cleanup API [\#6239](https://github.com/kokkos/kokkos/pull/6239)
- Add converting assignment to `DualView`:  [\#6474](https://github.com/kokkos/kokkos/pull/6474)


### Build System Changes

- Export `Kokkos_CXX_COMPILER_VERSION` [\#6282](https://github.com/kokkos/kokkos/pull/6282)
- Disable default oneDPL support in Trilinos [\#6342](https://github.com/kokkos/kokkos/pull/6342)

### Incompatibilities (i.e. breaking changes)
 - Ensure that `Kokkos::complex` only gets instantiated for cv-unqualified floating-point types  [\#6251](https://github.com/kokkos/kokkos/pull/6251)
 - Removed (deprecated-3) support for volatile join operators in reductions [\#6385](https://github.com/kokkos/kokkos/pull/6385)
 - Enforce `ViewCtorArgs` restrictions for `create_mirror_view` [\#6304](https://github.com/kokkos/kokkos/pull/6304)
 - SIMD types for ARM NEON are not autodetected anymore but need `Kokkos_ARCH_ARM_NEON` or `Kokkos_ARCH_NATIVE=ON` [\#6394](https://github.com/kokkos/kokkos/pull/6394)
 - Remove `#include <iostream>` from headers where possible [\#6482](https://github.com/kokkos/kokkos/pull/6482)

### Deprecations
- Deprecated `Kokkos::vector` [\#6252](https://github.com/kokkos/kokkos/pull/6252)
- All host allocation mechanisms except for `STD_MALLOC` have been deprecated [\#6341](https://github.com/kokkos/kokkos/pull/6341)

### Bug Fixes
 - Missing memory fence in `RandomPool::free_state` functions [\#6290](https://github.com/kokkos/kokkos/pull/6290)
 - Fix for corner case in `Kokkos::Experimental::is_partitioned` algorithm [\#6257](https://github.com/kokkos/kokkos/pull/6257)
 - Fix initialization of scratch lock variables in the `Cuda` backend [\#6433](https://github.com/kokkos/kokkos/pull/6433)
 - Fixes for `Kokkos::Array` [\#6372](https://github.com/kokkos/kokkos/pull/6372)
 - Fixed symlink configure issue for Windows [\#6241](https://github.com/kokkos/kokkos/pull/6241)
 - OpenMPTarget init-join fix [\#6444](https://github.com/kokkos/kokkos/pull/6444)
 - Fix atomic operations bug for Min and Max [\#6435](https://github.com/kokkos/kokkos/pull/6435)
 - Fix implementation for `cyl_bessel_i0` [\#6484](https://github.com/kokkos/kokkos/pull/6484)
 - Fix various NVCC warnings in `BinSort`, `Array`, and bit manipulation function templates [\#6483](https://github.com/kokkos/kokkos/pull/6483)

## [4.1.00](https://github.com/kokkos/kokkos/tree/4.1.00) (2023-06-16)
[Full Changelog](https://github.com/kokkos/kokkos/compare/4.0.01...4.1.00)

### Features:
* Add `<Kokkos_BitManipulation.hpp>` header [\#4577](https://github.com/kokkos/kokkos/pull/4577) [\#5907](https://github.com/kokkos/kokkos/pull/5907) [\#5967](https://github.com/kokkos/kokkos/pull/5967) [\#6101](https://github.com/kokkos/kokkos/pull/6101)
* Add `UnorderedMapInsertOpTypes` [\#5877](https://github.com/kokkos/kokkos/pull/5877) and documentation [\#350](https://github.com/kokkos/kokkos-core-wiki/pull/350)
* Add multiple reducers support for team-level parallel reduce [\#5727](https://github.com/kokkos/kokkos/pull/5727)

### Backend and Architecture Enhancements:

#### CUDA:

* Allow NVCC 12 to compile using C++20 flag [\#5977](https://github.com/kokkos/kokkos/pull/5977)
* Remove ability to disable CMake option `Kokkos_ENABLE_CUDA_LAMBDA` and unconditionally enable CUDA extended lambda support. [\#5964](https://github.com/kokkos/kokkos/pull/5964)
* Drop unnecessary fences around the memory allocation when using `CudaUVMSpace` in views [\#6008](https://github.com/kokkos/kokkos/pull/6008)

#### HIP:
* Improve performance for `parallel_reduce`. Use different parameters for `LightWeight` kernels [\#6029](https://github.com/kokkos/kokkos/pull/6029) and [\#6160](https://github.com/kokkos/kokkos/pull/6160)

#### SYCL:
* Only pass one wrapper object in SYCL reductions [\#6047](https://github.com/kokkos/kokkos/pull/6047)
* Improve and simplify parallel_scan implementation [\#6064](https://github.com/kokkos/kokkos/pull/6064)
* Remove workaround for submit_barrier not being enqueued properly [\#5504](https://github.com/kokkos/kokkos/pull/5504)
* Fix guards for using scratch space with SYCL [\#6003](https://github.com/kokkos/kokkos/pull/6003)
* Fix compiling SYCL with KOKKOS_IMPL_DO_NOT_USE_PRINTF_USAGE [\#6219](https://github.com/kokkos/kokkos/pull/6219)

#### OpenMPTarget:
* Improve hierarchical parallelism for Intel architectures [\#6043](https://github.com/kokkos/kokkos/pull/6043)
* Enable Cray compiler for the OpenMPTarget backend. [\#5889](https://github.com/kokkos/kokkos/pull/5889)

#### HPX:
* Update HPX backend to use HPX's sender/receiver functionality [\#5628](https://github.com/kokkos/kokkos/pull/5628)
* Increase minimum required HPX version to 1.8.0 [\#6132](https://github.com/kokkos/kokkos/pull/6132)
* Implement HPX::in_parallel [\#6143](https://github.com/kokkos/kokkos/pull/6143)

### General Enhancements
* Export CMake `Kokkos_{CUDA,HIP}_ARCHITECTURES` variables [\#5919](https://github.com/kokkos/kokkos/pull/5919) [\#5925](https://github.com/kokkos/kokkos/pull/5925)
* Add `Kokkos::Profiling::ScopedRegion` [\#5959](https://github.com/kokkos/kokkos/pull/5959) [\#5972](https://github.com/kokkos/kokkos/pull/5972)
* Add support for `View::rank[_dynamic]()`[\#5870](https://github.com/kokkos/kokkos/pull/5870)
* Detect incompatible relocatable device code mode to prevent ODR violations [\#5991](https://github.com/kokkos/kokkos/pull/5991)
* Add (experimental) support for 32-bit Darwin and PPC [\#5916](https://github.com/kokkos/kokkos/pull/5916)
* Add missing half and bhalf specialization of the infinity numeric trait [\#6055](https://github.com/kokkos/kokkos/pull/6055)
* Add `is_dual_view` trait and align further with regular view [\#6120](https://github.com/kokkos/kokkos/pull/6120)
* Allow templated functors in parallel_for, parallel_reduce and parallel_scan [\#5976](https://github.com/kokkos/kokkos/pull/5976)
* Define KOKKOS_COMPILER_INTEL_LLVM and only define at most one KOKKOS_COMPILER* macro [\#5906](https://github.com/kokkos/kokkos/pull/5906)
* Allow linking against build tree [\#6078](https://github.com/kokkos/kokkos/pull/6078)
* Allow passing a temporary std::vector to partition_space [\#6167](https://github.com/kokkos/kokkos/pull/6167)
* `Kokkos` can be used as an external dependency in `Trilinos` [\#6142](https://github.com/kokkos/kokkos/pull/6142), [\#6157](https://github.com/kokkos/kokkos/pull/6157) [\#6163](https://github.com/kokkos/kokkos/pull/6163)
* Left align demangled stacktrace output [\#6191](https://github.com/kokkos/kokkos/pull/6191)
* Improve OpenMP affinity warning to include MPI concerns [\#6185](https://github.com/kokkos/kokkos/pull/6185)

### Build System Changes
* Drop `Kokkos_ENABLE_LAUNCH_COMPILER` option which had no effect [\#6148](https://github.com/kokkos/kokkos/pull/6148)
* Export variables for relevant Kokkos options with cmake[\#6142](https://github.com/kokkos/kokkos/pull/6142)

### Incompatibilities (i.e. breaking changes)
* Desul atomics always enabled [\#5801](https://github.com/kokkos/kokkos/pull/5801)
* Drop `KOKKOS_ENABLE_CUDA_ASM*` and `KOKKOS_ENABLE_*_ATOMICS` macros [\#5940](https://github.com/kokkos/kokkos/pull/5940)
* Drop `KOKKOS_ENABLE_RFO_PREFETCH` macro [\#5944](https://github.com/kokkos/kokkos/pull/5944)
* Deprecate `Kokkos_ENABLE_CUDA_LAMBDA` configuration option and force it to `ON` [\#5964](https://github.com/kokkos/kokkos/pull/5964)
* Remove TriBITS Kokkos subpackages [\#6104](https://github.com/kokkos/kokkos/pull/6104)
* Cuda: Remove unused attach_texture_object [\#6129](https://github.com/kokkos/kokkos/pull/6129)
* Drop Kokkos_ENABLE_PROFILING_LOAD_PRINT configuration option [\#6150](https://github.com/kokkos/kokkos/pull/6150)
* Drop pointless Kokkos{Algorithms,Containers}_config.h files [\#6108](https://github.com/kokkos/kokkos/pull/6108)

### Deprecations
* Deprecate `BinSort`, `BinOp1D`, and `BinOp3D` default constructors [\#6131](https://github.com/kokkos/kokkos/pull/6131)

### Bug Fixes
* Fix `SYCLTeamMember` to take arguments for scratch sizes as `std::size_t` [\#5981](https://github.com/kokkos/kokkos/pull/5981)
* Fix Kokkos_SIMD with AVX2 on 64-bit architectures [\#6075](https://github.com/kokkos/kokkos/pull/6075)
* Fix an incorrectly returning size for SIMD uint64_t in AVX2 [\#6004](https://github.com/kokkos/kokkos/pull/6004)
* Fix missing avx512 header file with gcc versions before 10 [\#6183](https://github.com/kokkos/kokkos/pull/6183)
* Fix incorrect results of `parallel_reduce` of types smaller than `int` on CUDA and HIP: [\#5745](https://github.com/kokkos/kokkos/pull/5745)
* CMake: update package compatibility mode when building within Trilinos [\#6012](https://github.com/kokkos/kokkos/pull/6012)
* Fix warnings generated from internal uses of `ALL_t` rather than `Kokkos::ALL_t` [\#6028](https://github.com/kokkos/kokkos/pull/6028)
* Fix bug in `hpcbind` script: check for correct Slurm variable [\#6116](https://github.com/kokkos/kokkos/pull/6116)
* KokkosTools: Don't call callbacks before backends are initialized [\#6114](https://github.com/kokkos/kokkos/pull/6114)
* Fix global fence in Kokkos::resize(DynRankView) [\#6184](https://github.com/kokkos/kokkos/pull/6184)
* Fix `BinSort` support for strided views [\#6081](https://github.com/kokkos/kokkos/pull/6184)
* Fix missing `is_*_view` traits in containers [\#6195](https://github.com/kokkos/kokkos/pull/6195)
* Fix broken OpenMP target on NVHPC [\#6171](https://github.com/kokkos/kokkos/pull/6171)
* Sorting an empty view should exit early and not fail [\#6130](https://github.com/kokkos/kokkos/pull/6130)

## [4.0.01](https://github.com/kokkos/kokkos/tree/4.0.01) (2023-04-14)
[Full Changelog](https://github.com/kokkos/kokkos/compare/4.0.00...4.0.01)

### Backend and Architecture Enhancements:

#### CUDA:

- Allow NVCC 12 to compile using C++20 flag [\#6020](https://github.com/kokkos/kokkos/pull/6020)
- Add CUDA Ada architecture support [\#6022](https://github.com/kokkos/kokkos/pull/6022)

#### HIP:

- Add support for AMDGPU target NAVI31 / RX 7900 XT(X): gfx1100 [\#6021](https://github.com/kokkos/kokkos/pull/6021)
- HIP: Fix warning from `std::memcpy` [\#6019](https://github.com/kokkos/kokkos/pull/6019)

#### SYCL:
- Fix `SYCLTeamMember` to take arguments for scratch sizes as `std::size_t` [\#5986](https://github.com/kokkos/kokkos/pull/5986)

### General Enhancements
- Fixup 4.0 change log [\#6023](https://github.com/kokkos/kokkos/pull/6023)

### Build System Changes
- Cherry-pick TriBITS update from Trilinos [\#6037](https://github.com/kokkos/kokkos/pull/6037)
- CMake: update package compatibility mode when building within Trilinos [\#6013](https://github.com/kokkos/kokkos/pull/6013)

### Bug Fixes
- Fix an incorrectly returning size for SIMD uint64_t in AVX2 [\#6011](https://github.com/kokkos/kokkos/pull/6011)
- Desul atomics: wrong value for `desul::Impl::numeric_limits_max<uint64_t>` [\#6018](https://github.com/kokkos/kokkos/pull/6018)
- Fix warning in some user code when using std::memcpy [\#6000](https://github.com/kokkos/kokkos/pull/6000)
- Fix excessive build times using Makefile.kokkos [\#6068](https://github.com/kokkos/kokkos/pull/6068)

## [4.0.0](https://github.com/kokkos/kokkos/tree/4.0.00) (2023-02-21)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.7.01...4.0.00)

### Features:
- Allow value types without default constructor in `Kokkos::View` with `Kokkos::WithoutInitializing` [\#5307](https://github.com/kokkos/kokkos/pull/5307)
- `parallel_scan` with `View` as result type. [\#5146](https://github.com/kokkos/kokkos/pull/5146)
- Introduced `SharedSpace`, an alias for a `MemorySpace` that is accessible by every `ExecutionSpace`. The memory is moved and then accessed locally. [\#5289](https://github.com/kokkos/kokkos/pull/5289)
- Introduced `SharedHostPinnedSpace`, an alias for a `MemorySpace` that is accessible by every `ExecutionSpace`. The memory is pinned to the host and accessed via zero-copy access. [\#5405](https://github.com/kokkos/kokkos/pull/5405)
- Add team- and thread-level `sort`, `sort_by_key` algorithms. [\#5317](https://github.com/kokkos/kokkos/pull/5317)
- Groundwork for `MDSpan` integration. [\#4973](https://github.com/kokkos/kokkos/pull/4973) and [\#5304](https://github.com/kokkos/kokkos/pull/5304)
- Introduced MD version of hierarchical parallelism: `TeamThreadMDRange`, `ThreadVectorMDRange` and `TeamVectorMDRange`. [\#5238](https://github.com/kokkos/kokkos/pull/5238)

### Backend and Architecture Enhancements:

#### CUDA:
- Allow CUDA PTX forward compatibility [\#3612](https://github.com/kokkos/kokkos/pull/3612) [\#5536](https://github.com/kokkos/kokkos/pull/5536) [\#5527](https://github.com/kokkos/kokkos/pull/5527)
- Add support for NVIDIA Hopper GPU architecture [\#5538](https://github.com/kokkos/kokkos/pull/5538)
- Don't rely on synchronization behavior of default stream in CUDA and HIP [\#5391](https://github.com/kokkos/kokkos/pull/5391)
- Improve CUDA cache config settings [\#5706](https://github.com/kokkos/kokkos/pull/5706)

#### HIP:
 - Move `HIP`, `HIPSpace`, `HIPHostPinnedSpace`, and `HIPManagedSpace` out of the `Experimental` namespace [\#5383](https://github.com/kokkos/kokkos/pull/5383)
 - Don't rely on synchronization behavior of default stream in CUDA and HIP [\#5391](https://github.com/kokkos/kokkos/pull/5391)
 - Export AMD architecture flag when using Trilinos [\#5528](https://github.com/kokkos/kokkos/pull/5528)
 - Fix linking error (see [OLCF issue](https://docs.olcf.ornl.gov/systems/crusher_quick_start_guide.html#olcfdev-1167-kokkos-build-failures-with-prgenv-amd)) when using `amdclang`: [\#5539](https://github.com/kokkos/kokkos/pull/5539)
 - Remove support for MI25 and added support for Navi 1030 [\#5522](https://github.com/kokkos/kokkos/pull/5522)
 - Fix race condition when using `HSA_XNACK=1`  [\#5755](https://github.com/kokkos/kokkos/pull/5755)
 - Add parameter to force using GlobalMemory launch mechanism. This can be used when encountering compiler bugs with ROCm 5.3 and 5.4  [\#5796](https://github.com/kokkos/kokkos/pull/5796)

#### SYCL:
- Delegate choice of workgroup size for `parallel_reduce` with `RangePolicy` to the compiler. [\#5227](https://github.com/kokkos/kokkos/pull/5227)
- SYCL `RangePolicy`: manually specify workgroup size through chunk size [\#4875](https://github.com/kokkos/kokkos/pull/4875)

#### OpenMPTarget:
- Select the right device [\#5492](https://github.com/kokkos/kokkos/pull/5492)

#### OpenMP:
 - Add `partition_space` [\#5105](https://github.com/kokkos/kokkos/pull/5105)

### General Enhancements
- Implement `OffsetView` constructor taking `pair`s and `ViewCtorProp` [\#5303](https://github.com/kokkos/kokkos/pull/5303)
- Promote math constants to `Kokkos::numbers` namespace [\#5434](https://github.com/kokkos/kokkos/pull/5434)
- Add overloads of `hypot` math function that take 3 arguments [\#5341](https://github.com/kokkos/kokkos/pull/5341)
- Add `fma` fused multiply-add math function [\#5428](https://github.com/kokkos/kokkos/pull/5428)
- Views using `MemoryTraits::Atomic` don't need `volatile` overloads for the value type anymore. [\#5455](https://github.com/kokkos/kokkos/pull/5455)
- Added `is_team_handle` trait [\#5375](https://github.com/kokkos/kokkos/pull/5375)
- Refactor desul atomics to support compiling CUDA with NVC++ [\#5431](https://github.com/kokkos/kokkos/pull/5431) [\#5497](https://github.com/kokkos/kokkos/pull/5497) [\#5498](https://github.com/kokkos/kokkos/pull/5498)
- Support finding `libquadmath` with native compiler support [\#5286](https://github.com/kokkos/kokkos/pull/5286)
- Add architecture flags for MSVC [\#5673](https://github.com/kokkos/kokkos/pull/5673)
- SIMD backend for ARM NEON [\#5829](https://github.com/kokkos/kokkos/pull/5829)

### Build System Changes
- Let CMake determine OpenMP flags. [\#4105](https://github.com/kokkos/kokkos/pull/4105)
- Update minimum compiler versions. [\#5323](https://github.com/kokkos/kokkos/pull/5323)
- Makefile and CMake support for C++23 [\#5283](https://github.com/kokkos/kokkos/pull/5283)
- Do not add `-cuda` to the link line with NVHPC compiler when the CUDA backend is not actually enabled [\#5485](https://github.com/kokkos/kokkos/pull/5485)
- Only add `-latomic` in generated GNU makefiles when OpenMPTarget backend is enabled [\#5501](https://github.com/kokkos/kokkos/pull/5501) [\#5537](https://github.com/kokkos/kokkos/pull/5537) (3.7 patch release candidate)
- `Kokkos_ENABLE_CUDA_LAMBDA` now `ON` by default with NVCC [\#5580](https://github.com/kokkos/kokkos/pull/5580)
- Fix enabling of relocatable device code when using CUDA as CMake language [\#5564](https://github.com/kokkos/kokkos/pull/5564)
- Fix cmake configuration with CUDA 12 [\#5691](https://github.com/kokkos/kokkos/pull/5691)

### Incompatibilities (i.e. breaking changes)
- ***Require C++17***  [\#5277](https://github.com/kokkos/kokkos/pull/5277)
- Turn setting `Kokkos_CXX_STANDARD` into an error [\#5293](https://github.com/kokkos/kokkos/pull/5293)
- Remove all deprecations in Kokkos 3 [\#5297](https://github.com/kokkos/kokkos/pull/5297)
- Remove `KOKKOS_COMPILER_CUDA_VERSION` [\#5430](https://github.com/kokkos/kokkos/pull/5430)
- Drop `reciprocal_overflow_threshold` numeric trait [\#5326](https://github.com/kokkos/kokkos/pull/5326)
- Move `reduction_identity` out of `<Kokkos_NumericTraits.hpp>` into a new `<Kokkos_ReductionIdentity.hpp>` header [\#5450](https://github.com/kokkos/kokkos/pull/5450)
- Reduction and scan routines will report an error if the `join()` operator they would use takes `volatile`-qualified parameters [\#5409](https://github.com/kokkos/kokkos/pull/5409)
- `ENABLE_CUDA_UVM` is dropped in favor of using `SharedSpace` as `MemorySpace` explicitly [\#5608](https://github.com/kokkos/kokkos/pull/5608)
- Remove Kokkos_ENABLE_CUDA_LDG_INTRINSIC option [\#5623](https://github.com/kokkos/kokkos/pull/5623)
- Don't rely on synchronization behavior of default stream in CUDA and HIP - this potentially will break unintended implicit synchronization with other libraries such as MPI [\#5391](https://github.com/kokkos/kokkos/pull/5391)
- Make ExecutionSpace::concurrency() a non-static member function [\#5655](https://github.com/kokkos/kokkos/pull/5655) and related PRs
- Remove code guarded by `KOKKOS_ENABLE_DEPRECATED_CODE_3`

### Deprecations
- Deprecate `CudaUVMSpace::available()` which always returned `true` [\#5614](https://github.com/kokkos/kokkos/pull/5614)
- Deprecate `volatile`-qualified members from `Kokkos::pair` and `Kokkos::complex` [\#5412](https://github.com/kokkos/kokkos/pull/5412)
- Deprecate `KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_*` macros [\#5824](https://github.com/kokkos/kokkos/pull/5824) (oversight in 3.6)

### Bug Fixes
- Avoid allocating memory for `UniqueToken` [\#5300](https://github.com/kokkos/kokkos/pull/5300)
- Fix `pragma ivdep` in `Kokkos_OpenMP_Parallel.hpp` [\#5356](https://github.com/kokkos/kokkos/pull/5356)
- Fix configuring with Threads support when rerunning CMake [\#5486](https://github.com/kokkos/kokkos/pull/5486)
- Fix View assignment between `LayoutLeft` and `LayoutRight` with static extents [\#5535](https://github.com/kokkos/kokkos/pull/5535) (3.7 patch release candidate)
- Add `fence()` calls to sorting routine overloads that don't take an execution space parameter [\#5389](https://github.com/kokkos/kokkos/pull/5389)
- `ClockTic` changed to 64 bit to fix overflow on Power [\#5577](https://github.com/kokkos/kokkos/pull/5577) (incl. in 3.7.01 patch release)
- Fix incorrect offset in CUDA and HIP `parallel_scan` for < 4 byte types [\#5555](https://github.com/kokkos/kokkos/pull/5555) (3.7 patch release candidate)
- Fix incorrect alignment behavior of scratch allocations in some corner cases (e.g. very small allocations) [\#5687](https://github.com/kokkos/kokkos/pull/5687) (3.7 patch release candidate)
- Add missing `ReductionIdentity<char>` specialization [\#5798](https://github.com/kokkos/kokkos/pull/5798)
- Don't install standard algorithms headers multiple times [\#5670](https://github.com/kokkos/kokkos/pull/5670)
- Fix max scratch size calculation for level 0 scratch in CUDA and HIP [\#5718](https://github.com/kokkos/kokkos/pull/5718)

## [3.7.02](https://github.com/kokkos/kokkos/tree/3.7.02) (2023-05-17)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.7.01...3.7.02)

### Backends and Archs Enhancements:
#### CUDA
- Add Hopper support and update nvcc_wrapper to work with CUDA-12 [\#5693](https://github.com/kokkos/kokkos/pull/5693)
### General Enhancements:
- sprintf -> snprintf [\#5787](https://github.com/kokkos/kokkos/pull/5787)
### Build System:
- Add error message when not using `hipcc` and when `CMAKE_CXX_STANDARD` is not set [\#5945](https://github.com/kokkos/kokkos/pull/5945)
### Bug Fixes:
- Fix Scratch allocation alignment issues [\#5692](https://github.com/kokkos/kokkos/pull/5692)
- Fix Intel Classic Compiler ICE [\#5710](https://github.com/kokkos/kokkos/pull/5710)
- Don't install std algorithm headers multiple times [\#5711](https://github.com/kokkos/kokkos/pull/5711)
- Fix static init order issue in InitalizationSettings [\#5721](https://github.com/kokkos/kokkos/pull/5721)
- Fix src/dst Properties in deep_copy(DynamicView,View) [\#5732](https://github.com/kokkos/kokkos/pull/5732)
- Fix build on Fedora Rawhide [\#5782](https://github.com/kokkos/kokkos/pull/5782)
- Finalize HIP lock arrays [\#5694](https://github.com/kokkos/kokkos/pull/5694)
- Fix CUDA lock arrays for current Desul [\#5812](https://github.com/kokkos/kokkos/pull/5812)
- Set the correct device/context in InterOp tests [\#5701](https://github.com/kokkos/kokkos/pull/5701)

## [3.7.01](https://github.com/kokkos/kokkos/tree/3.7.01) (2022-12-01)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.7.00...3.7.01)

### Bug Fixes:
- Add fences to all sorting routines not taking an execution space instance argument [\#5547](https://github.com/kokkos/kokkos/pull/5547)
- Fix repeated `team_reduce` without barrier [\#5552](https://github.com/kokkos/kokkos/pull/5552)
- Fix memory spaces in `create_mirror_view` overloads using `view_alloc` [\#5521](https://github.com/kokkos/kokkos/pull/5521)
- Allow `as_view_of_rank_n()` to be overloaded for "special" scalar types [\#5553](https://github.com/kokkos/kokkos/pull/5553)
- Fix warning calling a `__host__` function from a `__host__ __device__` from `View:: as_view_of_rank_n` [\#5591](https://github.com/kokkos/kokkos/pull/5591)
- OpenMPTarget: adding implementation to set device id. [\#5557](https://github.com/kokkos/kokkos/pull/5557)
- Use `Kokkos::atomic_load` to Correct Race Condition Giving Rise to Seg Faulting Error in OpenMP tests [\#5559](https://github.com/kokkos/kokkos/pull/5559)
- cmake: define `KOKKOS_ARCH_A64FX` [\#5561](https://github.com/kokkos/kokkos/pull/5561)
- Only link against libatomic in gnu-make OpenMPTarget build [\#5565](https://github.com/kokkos/kokkos/pull/5565)
- Fix static extents assignment for LayoutLeft/LayoutRight assignment [\#5566](https://github.com/kokkos/kokkos/pull/5566)
- Do not add -cuda to the link line with NVHPC compiler when the CUDA backend is not actually enabled [\#5569](https://github.com/kokkos/kokkos/pull/5569)
- Export the flags in `KOKKOS_AMDGPU_OPTIONS` when using Trilinos [\#5571](https://github.com/kokkos/kokkos/pull/5571)
- Add support for detecting MPI local rank with MPICH and PMI [\#5570](https://github.com/kokkos/kokkos/pull/5570) [\#5582](https://github.com/kokkos/kokkos/pull/5582)
- Remove listing of undefined TPL dependencies [\#5573](https://github.com/kokkos/kokkos/pull/5573)
- ClockTic changed to 64 bit to fix overflow on Power [\#5592](https://github.com/kokkos/kokkos/pull/5592)
- Fix incorrect offset in CUDA and HIP parallel scan for < 4 byte types [\#5607](https://github.com/kokkos/kokkos/pull/5607)
- Fix initialization of Cuda lock arrays [\#5622](https://github.com/kokkos/kokkos/pull/5622)

## [3.7.00](https://github.com/kokkos/kokkos/tree/3.7.00) (2022-08-22)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.6.01...3.7.00)

### Features:
- Use non-volatile `join()` member functions and `operator+=` in `parallel_reduce/scan` [\#4931](https://github.com/kokkos/kokkos/pull/4931) [\#4954](https://github.com/kokkos/kokkos/pull/4954) [\#4951](https://github.com/kokkos/kokkos/pull/4951)
- Add `SIMD` sub package (requires C++17) [\#5016](https://github.com/kokkos/kokkos/pull/5016)
- Add `is_finalized()` [\#5247](https://github.com/kokkos/kokkos/pull/5247)
- Promote mathematical functions from `namespace Kokkos::Experimental` to `namespace Kokkos` [\#4791](https://github.com/kokkos/kokkos/pull/4791)
- Promote `min`, `max`, `clamp`, `minmax` functions from `namespace Kokkos::Experimental` to `namespace Kokkos` [\#5170](https://github.com/kokkos/kokkos/pull/5170)
- Add `round`, `logb`, `nextafter`, `copysign`, and `signbit` math functions [\#4768](https://github.com/kokkos/kokkos/pull/4768)
- Add `HIPManagedSpace`, similar to `CudaUVMSpace` [\#5112](https://github.com/kokkos/kokkos/pull/5112)
- Accept view construction allocation properties in `create_mirror[_view,_view_and_copy]` and `resize/realloc` [\#5125](https://github.com/kokkos/kokkos/pull/5125) [\#5095](https://github.com/kokkos/kokkos/pull/5095) [\#5035](https://github.com/kokkos/kokkos/pull/5035) [\#4805](https://github.com/kokkos/kokkos/pull/4805) [\#4844](https://github.com/kokkos/kokkos/pull/4844)
- Allow `MemorySpace::allocate()` to be called with execution space [\#4826](https://github.com/kokkos/kokkos/pull/4826)
- Experimental: Compile time view subscriber [\#4197](https://github.com/kokkos/kokkos/pull/4197)

### Backends and Archs Enhancements:
- Add support for Sapphire Rapids Intel architecture [\#5015](https://github.com/kokkos/kokkos/pull/5015)
- Add support for ICX, SKL and ICL Intel architectures [\#5013](https://github.com/kokkos/kokkos/pull/5013) [\#4929](https://github.com/kokkos/kokkos/pull/4929)
- Add arch flags for Intel GPU Ponte Vecchio [\#4932](https://github.com/kokkos/kokkos/pull/4932)
- SYCL: require GPU if GPU architecture was set at configuration time (i.e. do not allow fallback to CPU device) [\#5264](https://github.com/kokkos/kokkos/pull/5264) [\#5222](https://github.com/kokkos/kokkos/pull/5222)
- SYCL: Add `SYCL::sycl_queue()` for interoperability [\#5241](https://github.com/kokkos/kokkos/pull/5241)
- SYCL: Loosen restriction for using built-in `sycl::group_broadcast` [\#4552](https://github.com/kokkos/kokkos/pull/4552)
- SYCL: preserve address space [\#4396](https://github.com/kokkos/kokkos/pull/4396)
- OpenMPTarget: Adding a workaound for team scan [\#5219](https://github.com/kokkos/kokkos/pull/5219)
- OpenMPTarget: Adding logic to skip the kernel launch if `league_size=0` [\#5067](https://github.com/kokkos/kokkos/pull/5067)
- OpenMPTarget: Make sure `Kokkos::abort()` causes abnormal program termination when called on the host-side [\#4808](https://github.com/kokkos/kokkos/pull/4808)
- HIP: Make HIPHostPinnedSpace coarse-grained [\#5152](https://github.com/kokkos/kokkos/pull/5152)
- Refactor OpenMP `parallel_for` implementation to use more native OpenMP constructs [\#4664](https://github.com/kokkos/kokkos/pull/4664)
- Add option to optimize for local CPU architecture `Kokkos_ARCH_NATIVE` [\#4930](https://github.com/kokkos/kokkos/pull/4930)


### Implemented enhancements
- Add command line argument/environment variable to print the configuration [\#5233](https://github.com/kokkos/kokkos/pull/5233)
- Improve error message in view memory access violations [\#4950](https://github.com/kokkos/kokkos/pull/4950)
- Remove unnecessary fences in View initialization [\#4823](https://github.com/kokkos/kokkos/pull/4823)
- Make `View::shmem_size()` device-callable [\#4936](https://github.com/kokkos/kokkos/pull/4936)
- Update numerics support for `__float128` [\#5081](https://github.com/kokkos/kokkos/pull/5081)
- Add `log10` overload for `Kokkos::complex` [\#5009](https://github.com/kokkos/kokkos/pull/5009)
- Add `[[nodiscard]]` to `ScopeGuard` [\#5224](https://github.com/kokkos/kokkos/pull/5224)
- Add structured binding support for `Kokkos::Array` [\#4962](https://github.com/kokkos/kokkos/pull/4962)
- Enable accessing `Kokkos::Array` elements in constant expressions [\#4916](https://github.com/kokkos/kokkos/pull/4916)
- Mark `as_view_of_rank_n` as KOKKOS_FUNCTION [\#5248](https://github.com/kokkos/kokkos/pull/5248)
- Cleanup/rework fence overloads [\#5148](https://github.com/kokkos/kokkos/pull/5148)
- Assert that `Layout` construction from extents is valid in functions taking integer extents [\#5209](https://github.com/kokkos/kokkos/pull/5209)
- Add `fill_random` overload that takes an execution space as first argument [\#5181](https://github.com/kokkos/kokkos/pull/5181)
- Avoid some unnecessary fences in `parallel_reduce/scan` [\#5154](https://github.com/kokkos/kokkos/pull/5154)
- Include `KOKKOS_ENABLE_LIBDL` in options when printing configuration [\#5086](https://github.com/kokkos/kokkos/pull/5086)
- DynRankView: make `layout()` return the same as a corresponding static View [\#5026](https://github.com/kokkos/kokkos/pull/5026)
- Use `_mm_malloc` for icpx [\#5012](https://github.com/kokkos/kokkos/pull/5012)
- Avoid forcing matching execution spaces in `BinSort` constructor and `sort()` [\#4919](https://github.com/kokkos/kokkos/pull/4919)
- Check number of bins in `BinSort` [\#4890](https://github.com/kokkos/kokkos/pull/4890)
- Improve performance in parallel STL-like algorithms [\#4887](https://github.com/kokkos/kokkos/pull/4887) [\#4886](https://github.com/kokkos/kokkos/pull/4886)
- Disable `memset` on A64FX and launch `parallel_for` instead (performance) [\#4884](https://github.com/kokkos/kokkos/pull/4884)
- Allow non-power-of-two team sizes for team reductions and scans [\#4809](https://github.com/kokkos/kokkos/pull/4809)

#### Harmonization of Kokkos execution environment initialization:
- Warn when unable to detect local MPI rank and user explicitly asked for it [\#5263](https://github.com/kokkos/kokkos/pull/5263)
- Refactor parsing of command line arguments and environment variables [\#5221](https://github.com/kokkos/kokkos/pull/5221)
- Refactor device selection at initialization [\#5211](https://github.com/kokkos/kokkos/pull/5211)
- Rename tools settings for consistency [\#5201](https://github.com/kokkos/kokkos/pull/5201)
- Print help only once [\#5128](https://github.com/kokkos/kokkos/pull/5128)
- Update precedence rule in initialization [\#5130](https://github.com/kokkos/kokkos/pull/5130)
- Warn instead of just ignoring user settings when kokkos-tools is disabled [\#5088](https://github.com/kokkos/kokkos/pull/5088)
- Drop numa args in threads backend initialization [\#5127](https://github.com/kokkos/kokkos/pull/5127)
- Warn users when a flag prefixed with -[-]kokkos is not recognized and do not remove it [\#5256](https://github.com/kokkos/kokkos/pull/5256)
- Give back to Core what belongs to Core (aka moving tune_internals option from Tools back to Core) [\#5202](https://github.com/kokkos/kokkos/pull/5202)

#### Build system updates:
- `nvcc_wrapper`: filter out -pedantic-errors from nvcc options [\#5235](https://github.com/kokkos/kokkos/pull/5235)
- `nvcc_wrapper`: add known nvcc option --source-in-ptx [\#5052](https://github.com/kokkos/kokkos/pull/5052)
- Link libdl as interface library [\#5179](https://github.com/kokkos/kokkos/pull/5179)
- Only show GPU architectures with enabled corresponding backend [\#5119](https://github.com/kokkos/kokkos/pull/5119)
- Enable optional external desul build [\#5021](https://github.com/kokkos/kokkos/pull/5021) [\#5132](https://github.com/kokkos/kokkos/pull/5132)
- Export `Kokkos_CXX_STANDARD` variable with CMake [\#5068](https://github.com/kokkos/kokkos/pull/5068)
- Suppress warnings with nvc++ [\#5031](https://github.com/kokkos/kokkos/pull/5031)
- Disallow multiple host architectures in CMake [\#4996](https://github.com/kokkos/kokkos/pull/4996)
- Do not include compiler warning flags in the compile option of the cmake target [\#4989](https://github.com/kokkos/kokkos/pull/4989)
- AOT flags for OpenMPTarget targeting Intel GPUs [\#4915](https://github.com/kokkos/kokkos/pull/4915)
- Repurpose `Kokkos_ARCH_INTEL_GEN` for SYCL to mean JIT to be conforming with OMPT [\#4894](https://github.com/kokkos/kokkos/pull/4894)
- Replace amdgpu-target with offload-arch [\#4874](https://github.com/kokkos/kokkos/pull/4874)
- Do not enable `kokkos_launch_compiler` when `CMAKE_CXX_COMPILER_LAUNCHER` is set [\#4870](https://github.com/kokkos/kokkos/pull/4870)
- Move CMake version check up [\#4797](https://github.com/kokkos/kokkos/pull/4797)

### Incompatibilities:
- Remove `KOKKOS_THREAD_LOCAL` [\#5064](https://github.com/kokkos/kokkos/pull/5064)
- Remove `KOKKOS_ENABLE_POSIX_MEMALIGN` [\#5011](https://github.com/kokkos/kokkos/pull/5011)
- Remove unused `KOKKOS_ENABLE_TM` [\#4995](https://github.com/kokkos/kokkos/pull/4995)
- Remove unused cmakedefine `KOKKOS_ENABLE_COMPILER_WARNINGS` [\#4883](https://github.com/kokkos/kokkos/pull/4883)
- Remove unused `KOKKOS_ENABLE_DUALVIEW_MODIFY_CHECK` [\#4882](https://github.com/kokkos/kokkos/pull/4882)
- Drop Instruction Set Architecture (ISA) macros [\#4981](https://github.com/kokkos/kokkos/pull/4981)
- Warn in `ScopeGuard` about illegal usage [\#5250](https://github.com/kokkos/kokkos/pull/5250)

### Deprecations:
- Guard against non-public header inclusion [\#5178](https://github.com/kokkos/kokkos/pull/5178)
- Raise deprecation warnings if non empty WorkTag class is used [\#5230](https://github.com/kokkos/kokkos/pull/5230)
- Deprecate `parallel_*` overloads taking the label as trailing argument [\#5141](https://github.com/kokkos/kokkos/pull/5141)
- Deprecate nested types in functional [\#5185](https://github.com/kokkos/kokkos/pull/5185)
- Deprecate `InitArguments` struct and replace it with `InitializationSettings` [\#5135](https://github.com/kokkos/kokkos/pull/5135)
- Deprecate `finalize_all()` [\#5134](https://github.com/kokkos/kokkos/pull/5134)
- Deprecate command line arguments (other than `--help`) that are not prefixed with `kokkos-*` [\#5120](https://github.com/kokkos/kokkos/pull/5120)
- Deprecate `--[kokkos-]numa` cmdline arg and `KOKKOS_NUMA` env var [\#5117](https://github.com/kokkos/kokkos/pull/5117)
- Deprecate `--[kokkos-]threads` command line argument in favor of `--[kokkos-]num-threads` [\#5111](https://github.com/kokkos/kokkos/pull/5111)
- Deprecate `Kokkos::is_reducer_type` [\#4957](https://github.com/kokkos/kokkos/pull/4957)
- Deprecate `OffsetView` constructors taking `index_list_type` [\#4810](https://github.com/kokkos/kokkos/pull/4810)
- Deprecate overloads of `Kokkos::sort` taking a parameter `bool always_use_kokkos_sort` [\#5382](https://github.com/kokkos/kokkos/issues/5382)
- Warn about `parallel_reduce` cases that call `join()` with volatile-qualified arguments [\#5215](https://github.com/kokkos/kokkos/pull/5215)

### Bug Fixes:
- CUDA Reductions: Fix data races reported by Nvidia `compute-sanitizer` [\#4855](https://github.com/kokkos/kokkos/pull/4855)
- Work around Intel compiler bug [\#5301](https://github.com/kokkos/kokkos/pull/5301)
- Avoid allocating memory for UniqueToken [\#5300](https://github.com/kokkos/kokkos/pull/5300)
- DynamicView: Properly resize mirror instances after construction [\#5276](https://github.com/kokkos/kokkos/pull/5276)
- Remove Kokkos::Rank limit of 6 ranks [\#5271](https://github.com/kokkos/kokkos/pull/5271)
- Do not forget to set last element to nullptr when removing a flag in `Kokkos::initialize` [\#5272](https://github.com/kokkos/kokkos/pull/5272)
- Fix CUDA+MSVC build issue [\#5261](https://github.com/kokkos/kokkos/pull/5261)
- Fix `DynamicView::resize_serial` [\#5220](https://github.com/kokkos/kokkos/pull/5220)
- Fix cmake default compiler flags for unknown compiler [\#5217](https://github.com/kokkos/kokkos/pull/5217)
- Fix `move_backward` [\#5191](https://github.com/kokkos/kokkos/pull/5191)
- Fixing issue 5196 - missing symbol with intel compiler [\#5207](https://github.com/kokkos/kokkos/pull/5207)
- Preserve `KOKKOS_INVALID_INDEX` in ViewDimension and ArrayLayout construction [\#5188](https://github.com/kokkos/kokkos/pull/5188)
- Finalize `deep_copy_space` early avoiding printing to `std::cerr` for Cuda [\#5151](https://github.com/kokkos/kokkos/pull/5151)
- Use correct policy in Threads MDRange `parallel_reduce` [\#5123](https://github.com/kokkos/kokkos/pull/5123)
- Fix building with NVCC as the CXX compiler while the CUDA backend is not enabled [\#5115](https://github.com/kokkos/kokkos/pull/5115)
- OpenMPTarget Index range fix for MDRange. [\#5089](https://github.com/kokkos/kokkos/pull/5089)
- Fix bug with CUDA's team reduction for empty ranges [\#5079](https://github.com/kokkos/kokkos/pull/5079)
- Fix using `ZeroMemset` for Serial [\#5077](https://github.com/kokkos/kokkos/pull/5077)
- Fix `Kokkos::Vector::push_back` for default execution space [\#5047](https://github.com/kokkos/kokkos/pull/5047)
- ScatterView: Fix ScatterMin/ScatterMax to use proper atomics [\#5045](https://github.com/kokkos/kokkos/pull/5045)
- Fix calling `ZeroMemset` in `deep_copy` [\#5040](https://github.com/kokkos/kokkos/pull/5040)
- Make View self-assignment not produce double-free [\#5024](https://github.com/kokkos/kokkos/pull/5024)
- Guard against unrecognized pragma with intel compilers [\#5019](https://github.com/kokkos/kokkos/pull/5019)
- Fix racing condition in `HIPParallelLaunch` [\#5008](https://github.com/kokkos/kokkos/pull/5008)
- KokkosP: Fix `device_id` in profiling [\#4997](https://github.com/kokkos/kokkos/pull/4997)
- Fix for `Kokkos::vector::insert` into empty vector with begin and end iterators [\#4988](https://github.com/kokkos/kokkos/pull/4988)
- Fix Core header files installation [\#4984](https://github.com/kokkos/kokkos/pull/4984)
- Fix bounds errors with `Kokkos::sort` [\#4980](https://github.com/kokkos/kokkos/pull/4980)
- Fixup let `RangePolicy::set_chunk_size` return a reference to self [\#4918](https://github.com/kokkos/kokkos/pull/4918)
- Fix allocating large Views [\#4907](https://github.com/kokkos/kokkos/pull/4907)
- Fix combined reductions with `Kokkos::View` [\#4896](https://github.com/kokkos/kokkos/pull/4896)
- Fixed `_CUDA_ARCH__` to `__CUDA_ARCH__` for CUDA LDG [\#4893](https://github.com/kokkos/kokkos/pull/4893)
- Fixup `View::access()` truncate parameter pack [\#4876](https://github.com/kokkos/kokkos/pull/4876)
- Fix `abort` with HIP backend for ROCm 5.0.2 and beyond [\#4873](https://github.com/kokkos/kokkos/pull/4873)
- Fix HIP version when printing the configuration [\#4872](https://github.com/kokkos/kokkos/pull/4872)
- Fix scratch lock array when using scratch level 1 [\#4871](https://github.com/kokkos/kokkos/pull/4871)
- Fix Makefile.kokkos to work with fujitsu compiler [\#4867](https://github.com/kokkos/kokkos/pull/4867)
- cmake: Correct link THREADS link option [\#4854](https://github.com/kokkos/kokkos/pull/4854)
- UniqueToken `impl_acquire` function should be device only [\#4819](https://github.com/kokkos/kokkos/pull/4819)
- Fix example calls to non existing static `print_configuration` [\#4806](https://github.com/kokkos/kokkos/pull/4806)
- Fix requests for large team scratch sizes [\#4728](https://github.com/kokkos/kokkos/pull/4728)


## [3.6.01](https://github.com/kokkos/kokkos/tree/3.6.01) (2022-05-23)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.6.00...3.6.01)

### Bug Fixes:
- Fix Threads: Fix serial resizing scratch space (3.6.01 cherry-pick) [\#5109](https://github.com/kokkos/kokkos/pull/5109)
- Fix ScatterMin/ScatterMax to use proper atomics (3.6.01 cherry-pick) [\#5046](https://github.com/kokkos/kokkos/pull/5046)
- Fix allocating large Views [\#4907](https://github.com/kokkos/kokkos/pull/4907)
- Fix bounds errors with Kokkos::sort [\#4980](https://github.com/kokkos/kokkos/pull/4980)
- Fix HIP version when printing the configuration [\#4872](https://github.com/kokkos/kokkos/pull/4872)
- Fixed `_CUDA_ARCH__` to `__CUDA_ARCH__` for CUDA LDG [\#4893](https://github.com/kokkos/kokkos/pull/4893)
- Fixed an incorrect struct initialization [\#5028](https://github.com/kokkos/kokkos/pull/5028)
- Fix racing condition in `HIPParallelLaunch` [\#5008](https://github.com/kokkos/kokkos/pull/5008)
- Avoid deprecation warnings with `OpenMPExec::validate_partition` [\#4982](https://github.com/kokkos/kokkos/pull/4982)
- Make View self-assignment not produce double-free [\#5024](https://github.com/kokkos/kokkos/pull/5024)


## [3.6.00](https://github.com/kokkos/kokkos/tree/3.6.00) (2022-02-18)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.5.00...3.6.00)

### Features:
- Add C++ standard algorithms [\#4315](https://github.com/kokkos/kokkos/pull/4315)
- Implement `fill_random` for `DynRankView` [\#4763](https://github.com/kokkos/kokkos/pull/4763)
- Add `bhalf_t` [\#4543](https://github.com/kokkos/kokkos/pull/4543) [\#4653](https://github.com/kokkos/kokkos/pull/4653)
- Add mathematical constants [\#4519](https://github.com/kokkos/kokkos/pull/4519)
- Allow `Kokkos::{create_mirror*,resize,realloc}` to be used with `WithoutInitializing` [\#4486](https://github.com/kokkos/kokkos/pull/4486) [\#4337](https://github.com/kokkos/kokkos/pull/4337)
- Implement `KOKKOS_IF_ON_{HOST,DEVICE}` macros [\#4660](https://github.com/kokkos/kokkos/pull/4660)
- Allow setting the CMake language for Kokkos [\#4323](https://github.com/kokkos/kokkos/pull/4323)

#### Perf bug fix
- Desul: Add ScopeCaller [\#4690](https://github.com/kokkos/kokkos/pull/4690)
- Enable Desul atomics by default when using Makefiles [\#4606](https://github.com/kokkos/kokkos/pull/4606)
- Unique token improvement [\#4741](https://github.com/kokkos/kokkos/pull/4741) [\#4748](https://github.com/kokkos/kokkos/pull/4748)

#### Other improvements:
- Add math function long double overload on the host side [\#4712](https://github.com/kokkos/kokkos/pull/4712)

### Deprecations:
- Array reductions with pointer return types [\#4756](https://github.com/kokkos/kokkos/pull/4756)
- Deprecate `partition_master`, `validate_partition` [\#4737](https://github.com/kokkos/kokkos/pull/4737)
- Deprecate `Kokkos_ENABLE_PTHREAD` in favor of `Kokkos_ENABLE_THREADS` [\#4619](https://github.com/kokkos/kokkos/pull/4619) ** pair with use std::threads **
- Deprecate `log2(unsigned) -> int` (removing in next release) [\#4595](https://github.com/kokkos/kokkos/pull/4595)
- Deprecate `Kokkos::Impl::is_view` [\#4592](https://github.com/kokkos/kokkos/pull/4592)
- Deprecate `KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_*` macros and the `ActiveExecutionMemorySpace` alias [\#4668](https://github.com/kokkos/kokkos/issues/4668)

### Backends and Archs Enhancements:

#### SYCL:
- Update required SYCL compiler version [\#4749](https://github.com/kokkos/kokkos/pull/4749)
- Cap vector size to kernel maximum for SYCL [\#4704](https://github.com/kokkos/kokkos/pull/4704)
- Improve check for compatibility of vector size and subgroup size in SYCL [\#4579](https://github.com/kokkos/kokkos/pull/4579)
- Provide `chunk_size` for SYCL [\#4635](https://github.com/kokkos/kokkos/pull/4635)
- Use host-pinned memory for SYCL kernel memory [\#4627](https://github.com/kokkos/kokkos/pull/4627)
- Use shuffle-based algorithm for scalar reduction [\#4608](https://github.com/kokkos/kokkos/pull/4608)
- Implement pool of USM IndirectKernelMemory [\#4596](https://github.com/kokkos/kokkos/pull/4596)
- Provide valid default team size for SYCL [\#4481](https://github.com/kokkos/kokkos/pull/4481)

#### CUDA:
- Add checks for shmem usage in `parallel_reduce` [\#4548](https://github.com/kokkos/kokkos/pull/4548)

#### HIP:
- Add support for fp16 in the HIP backend [\#4688](https://github.com/kokkos/kokkos/pull/4688)
- Disable multiple kernel instantiations when using HIP (configure with `-DKokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS=ON` to use) [\#4644](https://github.com/kokkos/kokkos/pull/4644)
- Fix HIP scratch use per instance [\#4439](https://github.com/kokkos/kokkos/pull/4439)
- Change allocation header to 256B alignment for AMD VEGA architecture [\#4753](https://github.com/kokkos/kokkos/pull/4753)
- Add generic `KOKKOS_ARCH_VEGA` macro [\#4782](https://github.com/kokkos/kokkos/pull/4782)
- Require ROCm 4.5 [\#4689](https://github.com/kokkos/kokkos/pull/4689)

### HPX:
- Adapt to HPX 1.7.0 which is now required [\#4241](https://github.com/kokkos/kokkos/pull/4241)

#### OpenMP:
- Fix thread deduction for OpenMP for `thread_count==0` [\#4541](https://github.com/kokkos/kokkos/pull/4541)

#### OpenMPTarget:
- Update memory space `size_type` to improve performance (`size_t -> unsigned`) [\#4779](https://github.com/kokkos/kokkos/pull/4779)

#### Other Improvements:
- Improve NVHPC support [\#4599](https://github.com/kokkos/kokkos/pull/4599)
- Add `Kokkos::Experimental::{min,max,minmax,clamp}` [\#4629](https://github.com/kokkos/kokkos/pull/4629) [\#4506](https://github.com/kokkos/kokkos/pull/4506)
- Use device type as template argument in Containers and Algorithms [\#4724](https://github.com/kokkos/kokkos/pull/4724) [\#4675](https://github.com/kokkos/kokkos/pull/4675)
- Implement `Kokkos::sort` with execution space [\#4490](https://github.com/kokkos/kokkos/pull/4490)
- `Kokkos::resize` always error out for mismatch in runtime rank [\#4681](https://github.com/kokkos/kokkos/pull/4681)
- Print current call stack when calling `Kokkos::abort()` from the host [\#4672](https://github.com/kokkos/kokkos/pull/4672) [\#4671](https://github.com/kokkos/kokkos/pull/4671)
- Detect mismatch of execution spaces in functors [\#4655](https://github.com/kokkos/kokkos/pull/4655)
- Improve view label access on host [\#4647](https://github.com/kokkos/kokkos/pull/4647)
- Error out for `const` scalar return type in reduction [\#4645](https://github.com/kokkos/kokkos/pull/4645)
- Don't allow calling `UnorderdMap::value_at` for a set [\#4639](https://github.com/kokkos/kokkos/pull/4639)
- Add `KOKKOS_COMPILER_NVHPC` macro, disable `quiet_NaN` and `signaling_NaN` [\#4586](https://github.com/kokkos/kokkos/pull/4586)
- Improve performance of `local_deep_copy` [\#4511](https://github.com/kokkos/kokkos/pull/4511)
- Improve performance when sorting integers [\#4464](https://github.com/kokkos/kokkos/pull/4464)
- Add missing numeric traits (`denorm_min`, `reciprocal_overflow_threshold`, `{quiet,silent}_NaN}`) and make them work on cv-qualified types [\#4466](https://github.com/kokkos/kokkos/pull/4466) [\#4415](https://github.com/kokkos/kokkos/pull/4415) [\#4473](https://github.com/kokkos/kokkos/pull/4473) [\#4443](https://github.com/kokkos/kokkos/pull/4443)

### Implemented enhancements BuildSystem
- Manually compute IntelLLVM compiler version for older CMake versions [\#4760](https://github.com/kokkos/kokkos/pull/4760)
- Add Xptxas without = to `nvcc_wrapper` [\#4646](https://github.com/kokkos/kokkos/pull/4646)
- Use external GoogleTest optionally [\#4563](https://github.com/kokkos/kokkos/pull/4563)
- Silent warnings about multiple optimization flags with `nvcc_wrapper` [\#4502](https://github.com/kokkos/kokkos/pull/4502)
- Use the same flags in Makefile.kokkos for POWER7/8/9 as for CMake [\#4483](https://github.com/kokkos/kokkos/pull/4483)
- Fix support for A64FX architecture [\#4745](https://github.com/kokkos/kokkos/pull/4745)

### Incompatibilities:
- Drop `KOKKOS_ARCH_HIP` macro when using generated GNU makefiles [\#4786](https://github.com/kokkos/kokkos/pull/4786)
- Remove gcc-toolchain auto add for clang in Makefile.kokkos [\#4762](https://github.com/kokkos/kokkos/pull/4762)

### Bug Fixes:
- Lock constant memory in Cuda/HIP kernel launch with a mutex (thread safety) [\#4525](https://github.com/kokkos/kokkos/pull/4525)
- Fix overflow for large requested scratch allocation [\#4551](https://github.com/kokkos/kokkos/pull/4551)
- Fix Windows build in mingw [\#4564](https://github.com/kokkos/kokkos/pull/4564)
- Fix `kokkos_launch_compiler`: escape `$` character [\#4769](https://github.com/kokkos/kokkos/pull/4769) [\#4703](https://github.com/kokkos/kokkos/pull/4703)
- Fix math functions with NVCC and GCC 5 as host compiler [\#4733](https://github.com/kokkos/kokkos/pull/4733)
- Fix shared build with Intel19 [\#4725](https://github.com/kokkos/kokkos/pull/4725)
- Do not install empty `desul/src/` directory [\#4714](https://github.com/kokkos/kokkos/pull/4714)
- Fix wrong `device_id` computation in `identifier_from_devid` (Profiling Interface) [\#4694](https://github.com/kokkos/kokkos/pull/4694)
- Fix a bug in CUDA scratch memory pool (abnormally high memory consumption) [\#4673](https://github.com/kokkos/kokkos/pull/4673)
- Remove eval of command args in `hpcbind` [\#4630](https://github.com/kokkos/kokkos/pull/4630)
- SYCL fix to run when no GPU is detected [\#4623](https://github.com/kokkos/kokkos/pull/4623)
- Fix `layout_strides::span` for rank-0 views [\#4605](https://github.com/kokkos/kokkos/pull/4605)
- Fix SYCL atomics for local memory [\#4585](https://github.com/kokkos/kokkos/pull/4585)
- Hotfix `mdrange_large_deep_copy` for SYCL [\#4581](https://github.com/kokkos/kokkos/pull/4581)
- Fix bug when sorting integer using the HIP backend [\#4570](https://github.com/kokkos/kokkos/pull/4570)
- Fix compilation error when using HIP with RDC [\#4553](https://github.com/kokkos/kokkos/pull/4553)
- `DynamicView`: Fix deallocation extent [\#4533](https://github.com/kokkos/kokkos/pull/4533)
- SYCL fix running parallel_reduce with TeamPolicy for large ranges [\#4532](https://github.com/kokkos/kokkos/pull/4532)
- Fix bash syntax error in `nvcc_wrapper` [\#4524](https://github.com/kokkos/kokkos/pull/4524)
- OpenMPTarget `team_policy` reduce fixes for `init/join` reductions [\#4521](https://github.com/kokkos/kokkos/pull/4521)
- Avoid hangs in the Threads backend [\#4499](https://github.com/kokkos/kokkos/pull/4499)
- OpenMPTarget fix reduction bug in `parallel_reduce` for `TeamPolicy` [\#4491](https://github.com/kokkos/kokkos/pull/4491)
- HIP fix scratch space per instance [\#4439](https://github.com/kokkos/kokkos/pull/4439)
- OpenMPTarget fix team scratch allocation [\#4431](https://github.com/kokkos/kokkos/pull/4431)


## [3.5.00](https://github.com/kokkos/kokkos/tree/3.5.00) (2021-10-19)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.4.01...3.5.00)

### Features:

- Add support for quad-precision math functions/traits [\#4098](https://github.com/kokkos/kokkos/pull/4098)
- Adding ExecutionSpace partitioning function [\#4096](https://github.com/kokkos/kokkos/pull/4096)
- Improve Python Interop Capabilities [\#4065](https://github.com/kokkos/kokkos/pull/4065)
- Add half_t Kokkos::rand specialization [\#3922](https://github.com/kokkos/kokkos/pull/3922)
- Add math special functions: erf, erfcx, expint1, Bessel functions, Hankel functions [\#3920](https://github.com/kokkos/kokkos/pull/3920)
- Add missing common mathematical functions [\#4043](https://github.com/kokkos/kokkos/pull/4043) [\#4036](https://github.com/kokkos/kokkos/pull/4036) [\#4034](https://github.com/kokkos/kokkos/pull/4034)
- Let the numeric traits be SFINAE-friendly [\#4038](https://github.com/kokkos/kokkos/pull/4038)
- Add Desul atomics - enabling memory-order and memory-scope parameters [\#3247](https://github.com/kokkos/kokkos/pull/3247)
- Add detection idiom from the C++ standard library extension version 2 [\#3980](https://github.com/kokkos/kokkos/pull/3980)
- Fence Profiling Support in all backends [\#3966](https://github.com/kokkos/kokkos/pull/3966) [\#4304](https://github.com/kokkos/kokkos/pull/4304) [\#4258](https://github.com/kokkos/kokkos/pull/4258) [\#4232](https://github.com/kokkos/kokkos/pull/4232)
- Significant SYCL enhancements (see below)

### Deprecations:

- Deprecate CUDA_SAFE_CALL and HIP_SAFE_CALL [\#4249](https://github.com/kokkos/kokkos/pull/4249)
- Deprecate Kokkos::Impl::Timer (Kokkos::Timer has been available for a long time) [\#4201](https://github.com/kokkos/kokkos/pull/4201)
- Deprecate Experimental::MasterLock [\#4094](https://github.com/kokkos/kokkos/pull/4094)
- Deprecate Kokkos_TaskPolicy.hpp (headers got reorganized, doesn't remove functionality) [\#4011](https://github.com/kokkos/kokkos/pull/4011)
- Deprecate backward compatibility features [\#3978](https://github.com/kokkos/kokkos/pull/3978)
- Update and deprecate is_space::host_memory/execution/mirror_space [\#3973](https://github.com/kokkos/kokkos/pull/3973)


### Backends and Archs Enhancements:

- Enabling constbitset constructors in kernels [\#4296](https://github.com/kokkos/kokkos/pull/4296)
- Use ZeroMemset in View constructor to improve performance [\#4226](https://github.com/kokkos/kokkos/pull/4226)
- Use memset in deep_copy [\#3944](https://github.com/kokkos/kokkos/pull/3944)
- Add missing fence() calls in resize(View) that effectively do deep_copy(resized, orig) [\#4212](https://github.com/kokkos/kokkos/pull/4212)
- Avoid allocations in resize and realloc [\#4207](https://github.com/kokkos/kokkos/pull/4207)
- StaticCsrGraph: use device type instead of execution space to construct views [\#3991](https://github.com/kokkos/kokkos/pull/3991)
- Consider std::sort when view is accessible from host [\#3929](https://github.com/kokkos/kokkos/pull/3929)
- Fix CPP20 warnings except for volatile [\#4312](https://github.com/kokkos/kokkos/pull/4312)

#### SYCL:
- Introduce SYCLHostUSMSpace [\#4268](https://github.com/kokkos/kokkos/pull/4268)
- Implement SYCL TeamPolicy for vector_size > 1 [\#4183](https://github.com/kokkos/kokkos/pull/4183)
- Enable 64bit ranges for SYCL [\#4211](https://github.com/kokkos/kokkos/pull/4211)
- Don't print SYCL device info in execution space intialization [\#4168](https://github.com/kokkos/kokkos/pull/4168)
- Improve SYCL MDRangePolicy performance [\#4161](https://github.com/kokkos/kokkos/pull/4161)
- Use sub_groups in SYCL parallel_scan [\#4147](https://github.com/kokkos/kokkos/pull/4147)
- Implement subgroup reduction for SYCL RangePolicy parallel_reduce [\#3940](https://github.com/kokkos/kokkos/pull/3940)
- Use DPC++ broadcast extension in SYCL team_broadcast [\#4103](https://github.com/kokkos/kokkos/pull/4103)
- Only fence in SYCL parallel_reduce for non-device-accessible result_ptr [\#4089](https://github.com/kokkos/kokkos/pull/4089)
- Improve fencing behavior in SYCL backend [\#4088](https://github.com/kokkos/kokkos/pull/4088)
- Fence all registered SYCL queues before deallocating memory [\#4086](https://github.com/kokkos/kokkos/pull/4086)
- Implement SYCL::print_configuration [\#3992](https://github.com/kokkos/kokkos/pull/3992)
- Reuse scratch memory in parallel_scan and TeamPolicy (decreases memory footprint) [\#3899](https://github.com/kokkos/kokkos/pull/3899) [\#3889](https://github.com/kokkos/kokkos/pull/3889)

#### CUDA:
- Cuda improve heuristic for blocksize [\#4271](https://github.com/kokkos/kokkos/pull/4271)
- Don't use [[deprecated]] for nvcc [\#4229](https://github.com/kokkos/kokkos/pull/4229)
- Improve error message for NVHPC as host compiler [\#4227](https://github.com/kokkos/kokkos/pull/4227)
- Update support for cuda reductions to work with types < 4bytes [\#4156](https://github.com/kokkos/kokkos/pull/4156)
- Fix incompatible team size deduction in rare cases parallel_reduce [\#4142](https://github.com/kokkos/kokkos/pull/4142)
- Remove UVM usage in DynamicView [\#4129](https://github.com/kokkos/kokkos/pull/4129)
- Remove dependency between core and containers [\#4114](https://github.com/kokkos/kokkos/pull/4114)
- Adding opt-in CudaMallocSync support when using CUDA version >= 11.2 [\#4026](https://github.com/kokkos/kokkos/pull/4026) [\#4233](https://github.com/kokkos/kokkos/pull/4233)
- Fix a potential race condition in the CUDA backend [\#3999](https://github.com/kokkos/kokkos/pull/3999)

#### HIP:
- Implement new blocksize deduction method for HIP Backend [\#3953](https://github.com/kokkos/kokkos/pull/3953)
- Add multiple LaunchMechanism [\#3820](https://github.com/kokkos/kokkos/pull/3820)
- Make HIP backend thread-safe [\#4170](https://github.com/kokkos/kokkos/pull/4170)

#### Serial:
- Refactor Serial backend and fix thread-safety issue [\#4053](https://github.com/kokkos/kokkos/pull/4053)

#### OpenMPTarget:
- OpenMPTarget: support array reductions in RangePolicy [\#4040](https://github.com/kokkos/kokkos/pull/4040)
- OpenMPTarget: add MDRange parallel_reduce [\#4032](https://github.com/kokkos/kokkos/pull/4032)
- OpenMPTarget: Fix bug in for the case of a reducer. [\#4044](https://github.com/kokkos/kokkos/pull/4044)
- OpenMPTarget: verify process fix [\#4041](https://github.com/kokkos/kokkos/pull/4041)

### Implemented enhancements BuildSystem

#### Important BuildSystem Updates:
- Use hipcc architecture autodetection when Kokkos_ARCH is not set [\#3941](https://github.com/kokkos/kokkos/pull/3941)
- Introduce Kokkos_ENABLE_DEPRECATION_WARNINGS and remove deprecated code with Kokkos_ENABLE_DEPRECATED_CODE_3 [\#4106](https://github.com/kokkos/kokkos/pull/4106) [\#3855](https://github.com/kokkos/kokkos/pull/3855)

#### Other Improvements:
- Add allow-unsupported-compiler flag to nvcc-wrapper [\#4298](https://github.com/kokkos/kokkos/pull/4298)
- nvcc_wrapper: fix errors in argument handling [\#3993](https://github.com/kokkos/kokkos/pull/3993)
- Adds support for -time=<file> and -time <file> in nvcc_wrapper [\#4015](https://github.com/kokkos/kokkos/pull/4015)
- nvcc_wrapper: suppress duplicates of GPU architecture and RDC flags [\#3968](https://github.com/kokkos/kokkos/pull/3968)
- Fix TMPDIR support in nvcc_wrapper [\#3792](https://github.com/kokkos/kokkos/pull/3792)
- NVHPC: update PGI compiler arch flags [\#4133](https://github.com/kokkos/kokkos/pull/4133)
- Replace PGI with NVHPC (works for both) [\#4196](https://github.com/kokkos/kokkos/pull/4196)
- Make sure that KOKKOS_CXX_HOST_COMPILER_ID is defined [\#4235](https://github.com/kokkos/kokkos/pull/4235)
- Add options to Makefile builds for deprecated code and warnings [\#4215](https://github.com/kokkos/kokkos/pull/4215)
- Use KOKKOS_CXX_HOST_COMPILER_ID for identifying CPU arch flags [\#4199](https://github.com/kokkos/kokkos/pull/4199)
- Added support for Cray Clang to Makefile.kokkos [\#4176](https://github.com/kokkos/kokkos/pull/4176)
- Add XLClang as compiler [\#4120](https://github.com/kokkos/kokkos/pull/4120)
- Keep quoted compiler flags when passing to Trilinos [\#3987](https://github.com/kokkos/kokkos/pull/3987)
- Add support for AMD Zen3 CPU architecture [\#3972](https://github.com/kokkos/kokkos/pull/3972)
- Rename IntelClang to IntelLLVM [\#3945](https://github.com/kokkos/kokkos/pull/3945)
- Add cppcoreguidelines-pro-type-cstyle-cast to clang-tidy [\#3522](https://github.com/kokkos/kokkos/pull/3522)
- Add sve bit size definition for A64FX [\#3947](https://github.com/kokkos/kokkos/pull/3947) [\#3946](https://github.com/kokkos/kokkos/pull/3946)
- Remove KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES [\#4150](https://github.com/kokkos/kokkos/pull/4150)

### Other Changes:

#### Tool Enhancements:

- Retrieve original value from a point in a MultidimensionalSparseTuningProblem [\#3977](https://github.com/kokkos/kokkos/pull/3977)
- Allow extension of built-in tuners with additional tuning axes [\#3961](https://github.com/kokkos/kokkos/pull/3961)
- Added a categorical tuner [\#3955](https://github.com/kokkos/kokkos/pull/3955)


#### Miscellaneous:

- hpcbind: Use double quotes around $@ when invoking user command [\#4284](https://github.com/kokkos/kokkos/pull/4284)
- Add file and line to error message [\#3985](https://github.com/kokkos/kokkos/pull/3985)
- Fix compiler warnings when compiling with nvc++ [\#4198](https://github.com/kokkos/kokkos/pull/4198)
- Add OpenMPTarget CI build on AMD GPUs [\#4055](https://github.com/kokkos/kokkos/pull/4055)
- CI: icpx is now part of intel container [\#4002](https://github.com/kokkos/kokkos/pull/4002)

### Incompatibilities:

- Remove pre CUDA 9 KOKKOS_IMPL_CUDA_* macros [\#4138](https://github.com/kokkos/kokkos/pull/4138)

### Bug Fixes:
- UnorderedMap::clear() should zero the size() [\#4130](https://github.com/kokkos/kokkos/pull/4130)
- Add memory fence for HostSharedPtr::cleanup() [\#4144](https://github.com/kokkos/kokkos/pull/4144)
- SYCL: Fix race conditions in TeamPolicy::parallel_reduce [\#4418](https://github.com/kokkos/kokkos/pull/4418)
- Adding missing memory fence to serial exec space fence. [\#4292](https://github.com/kokkos/kokkos/pull/4292)
- Fix using external SYCL queues in tests [\#4291](https://github.com/kokkos/kokkos/pull/4291)
- Fix digits10 bug [\#4281](https://github.com/kokkos/kokkos/pull/4281)
- Fixes constexpr errors with frounding-math on gcc < 10. [\#4278](https://github.com/kokkos/kokkos/pull/4278)
- Fix compiler flags for PGI/NVHPC [\#4264](https://github.com/kokkos/kokkos/pull/4264)
- Fix Zen2/3 also implying Zen Arch with Makefiles [\#4260](https://github.com/kokkos/kokkos/pull/4260)
- Kokkos_Cuda.hpp: Fix shadow warning with cuda/11.0 [\#4252](https://github.com/kokkos/kokkos/pull/4252)
- Fix issue w/ static initialization of function attributes [\#4242](https://github.com/kokkos/kokkos/pull/4242)
- Disable long double hypot test on Power systems [\#4221](https://github.com/kokkos/kokkos/pull/4221)
- Fix false sharing in random pool [\#4218](https://github.com/kokkos/kokkos/pull/4218)
- Fix a missing memory_fence for debug shared alloc code [\#4216](https://github.com/kokkos/kokkos/pull/4216)
- Fix two xl issues [\#4179](https://github.com/kokkos/kokkos/pull/4179)
- Makefile.kokkos: fix (standard_in) 1: syntax error [\#4173](https://github.com/kokkos/kokkos/pull/4173)
- Fixes for query_device example [\#4172](https://github.com/kokkos/kokkos/pull/4172)
- Fix a bug when using HIP atomic with Kokkos::Complex [\#4159](https://github.com/kokkos/kokkos/pull/4159)
- Fix mistaken logic in pthread creation [\#4157](https://github.com/kokkos/kokkos/pull/4157)
- Define KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION when requesting Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=ON [\#4107](https://github.com/kokkos/kokkos/pull/4107)
- Fix compilation with latest MSVC version [\#4102](https://github.com/kokkos/kokkos/pull/4102)
- Fix incorrect macro definitions when compiling with Intel compiler on Windows [\#4087](https://github.com/kokkos/kokkos/pull/4087)
- Fixup global buffer overflow in hand rolled string manipulation [\#4070](https://github.com/kokkos/kokkos/pull/4070)
- Fixup heap buffer overflow in cmd line args parsing unit tests [\#4069](https://github.com/kokkos/kokkos/pull/4069)
- Only add quotes in compiler flags for Trilinos if necessary [\#4067](https://github.com/kokkos/kokkos/pull/4067)
- Fixed invocation of tools init callbacks [\#4061](https://github.com/kokkos/kokkos/pull/4061)
- Work around SYCL JIT compiler issues with static variables [\#4013](https://github.com/kokkos/kokkos/pull/4013)
- Fix TestDetectionIdiom.cpp test inclusion for Trilinos/TriBITS [\#4010](https://github.com/kokkos/kokkos/pull/4010)
- Fixup allocation headers with OpenMPTarget backend [\#4003](https://github.com/kokkos/kokkos/pull/4003)
- Add missing specialization for OMPT to Kokkos Random [\#3967](https://github.com/kokkos/kokkos/pull/3967)
- Disable hypot long double test on power arches [\#3962](https://github.com/kokkos/kokkos/pull/3962)
- Use different EBO workaround for MSVC (rebased) [\#3924](https://github.com/kokkos/kokkos/pull/3924)
- Fix SYCL Kokkos::Profiling::(de)allocateData calls [\#3928](https://github.com/kokkos/kokkos/pull/3928)

## [3.4.01](https://github.com/kokkos/kokkos/tree/3.4.01) (2021-05-19)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.4.00...3.4.01)

**Bug Fixes:**
- Windows: Remove atomic_compare_exchange_strong overload conflicts with Windows [\#4024](https://github.com/kokkos/kokkos/pull/4024)
- OpenMPTarget: Fixup allocation headers with OpenMPTarget backend [\#4020](https://github.com/kokkos/kokkos/pull/4020)
- OpenMPTarget: Add missing specailization for OMPT to Kokkos Random [\#4022](https://github.com/kokkos/kokkos/pull/4022)
- AMD: Add support for AMD Zen3 CPU architecture [\#4021](https://github.com/kokkos/kokkos/pull/4021)
- SYCL: Implement SYCL::print_configuration [\#4012](https://github.com/kokkos/kokkos/pull/4012)
- Containers: staticcsrgraph: use device type instead of execution space to construct views [\#3998](https://github.com/kokkos/kokkos/pull/3998)
- nvcc_wrapper: fix errors in argument handling, suppress duplicates of GPU architecture and RDC flags [\#4006](https://github.com/kokkos/kokkos/pull/4006)
- CI: Add icpx testing to intel container [\#4004](https://github.com/kokkos/kokkos/pull/4004)
- CMake/TRIBITS: Keep quoted compiler flags when passing to Trilinos [\#4007](https://github.com/kokkos/kokkos/pull/4007)
- CMake: Rename IntelClang to IntelLLVM [\#3945](https://github.com/kokkos/kokkos/pull/3945)

## [3.4.00](https://github.com/kokkos/kokkos/tree/3.4.00) (2021-04-25)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.3.01...3.4.00)

**Highlights:**
- SYCL Backend Almost Feature Complete
- OpenMPTarget Backend Almost Feature Complete
- Performance Improvements for HIP backend
- Require CMake 3.16 or newer
- Tool Callback Interface Enhancements
- cmath wrapper functions available now in Kokkos::Experimental

**Features:**
- Implement parallel_scan with ThreadVectorRange and Reducer [\#3861](https://github.com/kokkos/kokkos/pull/3861)
- Implement SYCL Random [\#3849](https://github.com/kokkos/kokkos/pull/3849)
- OpenMPTarget: Adding Implementation for nested reducers [\#3845](https://github.com/kokkos/kokkos/pull/3845)
- Implement UniqueToken for SYCL [\#3833](https://github.com/kokkos/kokkos/pull/3833)
- OpenMPTarget: UniqueToken::Global implementation [\#3823](https://github.com/kokkos/kokkos/pull/3823)
- DualView sync's on ExecutionSpaces [\#3822](https://github.com/kokkos/kokkos/pull/3822)
- SYCL outer TeamPolicy parallel_reduce [\#3818](https://github.com/kokkos/kokkos/pull/3818)
- SYCL TeamPolicy::team_scan [\#3815](https://github.com/kokkos/kokkos/pull/3815)
- SYCL MDRangePolicy parallel_reduce [\#3801](https://github.com/kokkos/kokkos/pull/3801)
- Enable use of execution space instances in ScatterView [\#3786](https://github.com/kokkos/kokkos/pull/3786)
- SYCL TeamPolicy nested parallel_reduce [\#3783](https://github.com/kokkos/kokkos/pull/3783)
- OpenMPTarget: MDRange with TagType for parallel_for [\#3781](https://github.com/kokkos/kokkos/pull/3781)
- Adding OpenMPTarget parallel_scan [\#3655](https://github.com/kokkos/kokkos/pull/3655)
- SYCL basic TeamPolicy [\#3654](https://github.com/kokkos/kokkos/pull/3654)
- OpenMPTarget: scratch memory implementation [\#3611](https://github.com/kokkos/kokkos/pull/3611)

**Implemented enhancements Backends and Archs:**
- SYCL choose a specific GPU [\#3918](https://github.com/kokkos/kokkos/pull/3918)
- [HIP] Lock access to scratch memory when using Teams [\#3916](https://github.com/kokkos/kokkos/pull/3916)
- [HIP] fix multithreaded access to get_next_driver [\#3908](https://github.com/kokkos/kokkos/pull/3908)
- Forward declare HIPHostPinnedSpace and SYCLSharedUSMSpace [\#3902](https://github.com/kokkos/kokkos/pull/3902)
- Let SYCL USMObjectMem use SharedAllocationRecord [\#3898](https://github.com/kokkos/kokkos/pull/3898)
- Implement clock_tic for SYCL [\#3893](https://github.com/kokkos/kokkos/pull/3893)
- Don't use a static variable in HIPInternal::scratch_space [\#3866](https://github.com/kokkos/kokkos/pull/3866)(https://github.com/kokkos/kokkos/pull/3866)
- Reuse memory for SYCL parallel_reduce [\#3873](https://github.com/kokkos/kokkos/pull/3873)
- Update SYCL compiler in CI [\#3826](https://github.com/kokkos/kokkos/pull/3826)
- Introduce HostSharedPtr to manage m_space_instance for Cuda/HIP/SYCL [\#3824](https://github.com/kokkos/kokkos/pull/3824)
- [HIP] Use shuffle for range reduction [\#3811](https://github.com/kokkos/kokkos/pull/3811)
- OpenMPTarget: Changes to the hierarchical parallelism [\#3808](https://github.com/kokkos/kokkos/pull/3808)
- Remove ExtendedReferenceWrapper for SYCL parallel_reduce [\#3802](https://github.com/kokkos/kokkos/pull/3802)
- Eliminate sycl_indirect_launch [\#3777](https://github.com/kokkos/kokkos/pull/3777)
- OpenMPTarget: scratch implementation for parallel_reduce [\#3776](https://github.com/kokkos/kokkos/pull/3776)
- Allow initializing SYCL execution space from sycl::queue and SYCL::impl_static_fence [\#3767](https://github.com/kokkos/kokkos/pull/3767)
- SYCL TeamPolicy scratch memory alternative [\#3763](https://github.com/kokkos/kokkos/pull/3763)
- Alternative implementation for SYCL TeamPolicy [\#3759](https://github.com/kokkos/kokkos/pull/3759)
- Unify handling of synchronous errors in SYCL [\#3754](https://github.com/kokkos/kokkos/pull/3754)
- core/Cuda: Half_t updates for cgsolve [\#3746](https://github.com/kokkos/kokkos/pull/3746)
- Unify HIPParallelLaunch structures [\#3733](https://github.com/kokkos/kokkos/pull/3733)
- Improve performance for SYCL parallel_reduce [\#3732](https://github.com/kokkos/kokkos/pull/3732)
- Use consistent types in Kokkos_OpenMPTarget_Parallel.hpp [\#3703](https://github.com/kokkos/kokkos/pull/3703)
- Implement non-blocking kernel launches for HIP backend [\#3697](https://github.com/kokkos/kokkos/pull/3697)
- Change SYCLInternal::m_queue std::unique_ptr -> std::optional [\#3677](https://github.com/kokkos/kokkos/pull/3677)
- Use alternative SYCL parallel_reduce implementation [\#3671](https://github.com/kokkos/kokkos/pull/3671)
- Use runtime values in KokkosExp_MDRangePolicy.hpp [\#3626](https://github.com/kokkos/kokkos/pull/3626)
- Clean up AnalyzePolicy [\#3564](https://github.com/kokkos/kokkos/pull/3564)
- Changes for indirect launch of SYCL parallel reduce [\#3511](https://github.com/kokkos/kokkos/pull/3511)

**Implemented enhancements BuildSystem:**
- Also require C++14 when building gtest [\#3912](https://github.com/kokkos/kokkos/pull/3912)
- Fix compiling SYCL with OpenMP [\#3874](https://github.com/kokkos/kokkos/pull/3874)
- Require C++17 for SYCL (at configuration time) [\#3869](https://github.com/kokkos/kokkos/pull/3869)
- Add COMPILE_DEFINITIONS argument to kokkos_create_imported_tpl [\#3862](https://github.com/kokkos/kokkos/pull/3862)
- Do not pass arch flags to the linker with no rdc [\#3846](https://github.com/kokkos/kokkos/pull/3846)
- Try compiling C++14 check with C++14 support and print error message [\#3843](https://github.com/kokkos/kokkos/pull/3843)
- Enable HIP with Cray Clang [\#3842](https://github.com/kokkos/kokkos/pull/3842)
- Add an option to disable header self containment tests [\#3834](https://github.com/kokkos/kokkos/pull/3834)
- CMake check for C++14 [\#3809](https://github.com/kokkos/kokkos/pull/3809)
- Prefer -std=* over --std=* [\#3779](https://github.com/kokkos/kokkos/pull/3779)
- Kokkos launch compiler updates [\#3778](https://github.com/kokkos/kokkos/pull/3778)
- Updated comments and enabled no-op for kokkos_launch_compiler [\#3774](https://github.com/kokkos/kokkos/pull/3774)
- Apple's Clang not correctly recognised [\#3772](https://github.com/kokkos/kokkos/pull/3772)
- kokkos_launch_compiler + CUDA auto-detect arch [\#3770](https://github.com/kokkos/kokkos/pull/3770)
- Add Spack test support for Kokkos [\#3753](https://github.com/kokkos/kokkos/pull/3753)
- Split SYCL tests for aot compilation [\#3741](https://github.com/kokkos/kokkos/pull/3741)
- Use consistent OpenMP flag for IntelClang [\#3735](https://github.com/kokkos/kokkos/pull/3735)
- Add support for -Wno-deprecated-gpu-targets [\#3722](https://github.com/kokkos/kokkos/pull/3722)
- Add configuration to target CUDA compute capability 8.6 [\#3713](https://github.com/kokkos/kokkos/pull/3713)
- Added VERSION and SOVERSION to KOKKOS_INTERNAL_ADD_LIBRARY [\#3706](https://github.com/kokkos/kokkos/pull/3706)
- Add fast-math to known NVCC flags [\#3699](https://github.com/kokkos/kokkos/pull/3699)
- Add MI-100 arch string [\#3698](https://github.com/kokkos/kokkos/pull/3698)
- Require CMake >=3.16 [\#3679](https://github.com/kokkos/kokkos/pull/3679)
- KokkosCI.cmake, KokkosCTest.cmake.in, CTestConfig.cmake.in + CI updates [\#2844](https://github.com/kokkos/kokkos/pull/2844)

**Implemented enhancements Tools:**
- Improve readability of the callback invocation in profiling [\#3860](https://github.com/kokkos/kokkos/pull/3860)
- V1.1 Tools Interface: incremental, action-based [\#3812](https://github.com/kokkos/kokkos/pull/3812)
- Enable launch latency simulations [\#3721](https://github.com/kokkos/kokkos/pull/3721)
- Added metadata callback to tools interface [\#3711](https://github.com/kokkos/kokkos/pull/3711)
- MDRange Tile Size Tuning [\#3688](https://github.com/kokkos/kokkos/pull/3688)
- Added support for command-line args for kokkos-tools [\#3627](https://github.com/kokkos/kokkos/pull/3627)
- Query max tile sizes for an MDRangePolicy, and set tile sizes on an existing policy [\#3481](https://github.com/kokkos/kokkos/pull/3481)

**Implemented enhancements Other:**
- Try detecting ndevices in get_gpu [\#3921](https://github.com/kokkos/kokkos/pull/3921)
- Use strcmp to compare names() [\#3909](https://github.com/kokkos/kokkos/pull/3909)
- Add execution space arguments for constructor overloads that might allocate a new underlying View [\#3904](https://github.com/kokkos/kokkos/pull/3904)
- Prefix labels in internal use of kokkos_malloc [\#3891](https://github.com/kokkos/kokkos/pull/3891)
- Prefix labels for internal uses of SharedAllocationRecord [\#3890](https://github.com/kokkos/kokkos/pull/3890)
- Add missing hypot math function [\#3880](https://github.com/kokkos/kokkos/pull/3880)
- Unify algorithm unit tests to avoid code duplication [\#3851](https://github.com/kokkos/kokkos/pull/3851)
- DualView.template view() better matches for Devices in UVMSpace cases [\#3857](https://github.com/kokkos/kokkos/pull/3857)
- More extensive disentangling of Policy Traits [\#3829](https://github.com/kokkos/kokkos/pull/3829)
- Replaced nanosleep and sched_yield with STL routines [\#3825](https://github.com/kokkos/kokkos/pull/3825)
- Constructing Atomic Subviews [\#3810](https://github.com/kokkos/kokkos/pull/3810)
- Metadata Declaration in Core [\#3729](https://github.com/kokkos/kokkos/pull/3729)
- Allow using tagged final functor in parallel_reduce [\#3714](https://github.com/kokkos/kokkos/pull/3714)
- Major duplicate code removal in SharedAllocationRecord specializations [\#3658](https://github.com/kokkos/kokkos/pull/3658)

**Fixed bugs:**
- Provide forward declarations in Kokkos_ViewLayoutTiled.hpp for XL [\#3911](https://github.com/kokkos/kokkos/pull/3911)
- Fixup absolute value of floating points in Kokkos complex [\#3882](https://github.com/kokkos/kokkos/pull/3882)
- Address intel 17 ICE [\#3881](https://github.com/kokkos/kokkos/pull/3881)
- Add missing pow(Kokkos::complex) overloads [\#3868](https://github.com/kokkos/kokkos/pull/3868)
- Fix bug {pow, log}(Kokkos::complex) [\#3866](https://github.com/kokkos/kokkos/pull/3866)(https://github.com/kokkos/kokkos/pull/3866)
- Cleanup writing to output streams in Cuda [\#3859](https://github.com/kokkos/kokkos/pull/3859)
- Fixup cache CUDA fallback execution space instance used by DualView::sync [\#3856](https://github.com/kokkos/kokkos/pull/3856)
- Fix cmake warning with pthread [\#3854](https://github.com/kokkos/kokkos/pull/3854)
- Fix typo FOUND_CUDA_{DRIVVER -> DRIVER} [\#3852](https://github.com/kokkos/kokkos/pull/3852)
- Fix bug in SYCL team_reduce [\#3848](https://github.com/kokkos/kokkos/pull/3848)
- Atrocious bug in MDRange tuning [\#3803](https://github.com/kokkos/kokkos/pull/3803)
- Fix compiling SYCL with Kokkos_ENABLE_TUNING=ON [\#3800](https://github.com/kokkos/kokkos/pull/3800)
- Fixed command line parsing bug [\#3797](https://github.com/kokkos/kokkos/pull/3797)
- Workaround race condition in SYCL parallel_reduce [\#3782](https://github.com/kokkos/kokkos/pull/3782)
- Fix Atomic{Min,Max} for Kepler30 [\#3780](https://github.com/kokkos/kokkos/pull/3780)
- Fix SYCL typo [\#3755](https://github.com/kokkos/kokkos/pull/3755)
- Fixed Kokkos_install_additional_files macro [\#3752](https://github.com/kokkos/kokkos/pull/3752)
- Fix a typo for Kokkos_ARCH_A64FX [\#3751](https://github.com/kokkos/kokkos/pull/3751)
- OpenMPTarget: fixes and workarounds to work with "Release" build type [\#3748](https://github.com/kokkos/kokkos/pull/3748)
- Fix parsing bug for number of devices command line argument [\#3724](https://github.com/kokkos/kokkos/pull/3724)
- Avoid more warnings with clang and C++20 [\#3719](https://github.com/kokkos/kokkos/pull/3719)
- Fix gcc-10.1 C++20 warnings [\#3718](https://github.com/kokkos/kokkos/pull/3718)
- Fix cuda cache config not being set correct [\#3712](https://github.com/kokkos/kokkos/pull/3712)
- Fix dualview deepcopy perftools [\#3701](https://github.com/kokkos/kokkos/pull/3701)
- use drand instead of frand in drand [\#3696](https://github.com/kokkos/kokkos/pull/3696)

**Incompatibilities:**
- Remove unimplemented member functions of SYCLDevice [\#3919](https://github.com/kokkos/kokkos/pull/3919)
- Replace cl::sycl [\#3896](https://github.com/kokkos/kokkos/pull/3896)
- Get rid of SYCL workaround in Kokkos_Complex.hpp [\#3884](https://github.com/kokkos/kokkos/pull/3884)
- Replace most uses of if_c [\#3883](https://github.com/kokkos/kokkos/pull/3883)
- Remove Impl::enable_if_type [\#3863](https://github.com/kokkos/kokkos/pull/3863)
- Remove HostBarrier test [\#3847](https://github.com/kokkos/kokkos/pull/3847)
- Avoid (void) interface [\#3836](https://github.com/kokkos/kokkos/pull/3836)
- Remove VerifyExecutionCanAccessMemorySpace [\#3813](https://github.com/kokkos/kokkos/pull/3813)
- Avoid duplicated code in ScratchMemorySpace [\#3793](https://github.com/kokkos/kokkos/pull/3793)
- Remove superfluous FunctorFinal specialization [\#3788](https://github.com/kokkos/kokkos/pull/3788)
- Rename cl::sycl -> sycl in Kokkos_MathematicalFunctions.hpp [\#3678](https://github.com/kokkos/kokkos/pull/3678)
- Remove integer_sequence backward compatibility implementation [\#3533](https://github.com/kokkos/kokkos/pull/3533)

**Enabled tests:**
- Fixup re-enable core performance tests [\#3903](https://github.com/kokkos/kokkos/pull/3903)
- Enable more SYCL tests [\#3900](https://github.com/kokkos/kokkos/pull/3900)
- Restrict MDRange Policy tests for Intel GPUs [\#3853](https://github.com/kokkos/kokkos/pull/3853)
- Disable death tests for rawhide [\#3844](https://github.com/kokkos/kokkos/pull/3844)
- OpenMPTarget: Block unit tests that do not pass with the nvidia compiler [\#3839](https://github.com/kokkos/kokkos/pull/3839)
- Enable Bitset container test for SYCL [\#3830](https://github.com/kokkos/kokkos/pull/3830)
- Enable some more SYCL tests [\#3744](https://github.com/kokkos/kokkos/pull/3744)
- Enable SYCL atomic tests [\#3742](https://github.com/kokkos/kokkos/pull/3742)
- Enable more SYCL perf_tests [\#3692](https://github.com/kokkos/kokkos/pull/3692)
- Enable examples for SYCL [\#3691](https://github.com/kokkos/kokkos/pull/3691)

## [3.3.01](https://github.com/kokkos/kokkos/tree/3.3.01) (2021-01-06)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.3.00...3.3.01)

**Bug Fixes:**
- Fix severe performance bug in DualView which added memcpys for sync and modify [\#3693](https://github.com/kokkos/kokkos/issues/#3693)
- Fix performance bug in CUDA backend, where the cuda Cache config was not set correct.

## [3.3.00](https://github.com/kokkos/kokkos/tree/3.3.00) (2020-12-16)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.2.01...3.3.00)

**Features:**
- Require C++14 as minimum C++ standard. C++17 and C++20 are supported too.
- HIP backend is nearly feature complete. Kokkos Dynamic Task Graphs are missing.
- Major update for OpenMPTarget: many capabilities now work. For details contact us.
- Added DPC++/SYCL backend: primary capabilites are working.
- Added Kokkos Graph API analogous to CUDA Graphs.
- Added parallel_scan support with TeamThreadRange [\#3536](https://github.com/kokkos/kokkos/pull/3536)
- Added Logical Memory Spaces [\#3546](https://github.com/kokkos/kokkos/pull/3546)
- Added initial half precision support [\#3439](https://github.com/kokkos/kokkos/pull/3439)
- Experimental feature: control cuda occupancy [\#3379](https://github.com/kokkos/kokkos/pull/3379)

**Implemented enhancements Backends and Archs:**
- Add a64fx and fujitsu Compiler support [\#3614](https://github.com/kokkos/kokkos/pull/3614)
- Adding support for AMD gfx908 archictecture [\#3375](https://github.com/kokkos/kokkos/pull/3375)
- SYCL parallel\_for MDRangePolicy [\#3583](https://github.com/kokkos/kokkos/pull/3583)
- SYCL add parallel\_scan [\#3577](https://github.com/kokkos/kokkos/pull/3577)
- SYCL custom reductions [\#3544](https://github.com/kokkos/kokkos/pull/3544)
- SYCL Enable container unit tests [\#3550](https://github.com/kokkos/kokkos/pull/3550)
- SYCL feature level 5 [\#3480](https://github.com/kokkos/kokkos/pull/3480)
- SYCL Feature level 4 (parallel\_for) [\#3474](https://github.com/kokkos/kokkos/pull/3474)
- SYCL feature level 3 [\#3451](https://github.com/kokkos/kokkos/pull/3451)
- SYCL feature level 2 [\#3447](https://github.com/kokkos/kokkos/pull/3447)
- OpenMPTarget: Hierarchial reduction for + operator on scalars [\#3504](https://github.com/kokkos/kokkos/pull/3504)
- OpenMPTarget hierarchical [\#3411](https://github.com/kokkos/kokkos/pull/3411)
- HIP Add Impl::atomic\_[store,load] [\#3440](https://github.com/kokkos/kokkos/pull/3440)
- HIP enable global lock arrays [\#3418](https://github.com/kokkos/kokkos/pull/3418)
- HIP Implement multiple occupancy paths for various HIP kernel launchers [\#3366](https://github.com/kokkos/kokkos/pull/3366)

**Implemented enhancements Policies:**
- MDRangePolicy: Let it be semiregular [\#3494](https://github.com/kokkos/kokkos/pull/3494)
- MDRangePolicy: Check narrowing conversion in construction [\#3527](https://github.com/kokkos/kokkos/pull/3527)
- MDRangePolicy: CombinedReducers support [\#3395](https://github.com/kokkos/kokkos/pull/3395)
- Kokkos Graph: Interface and Default Implementation [\#3362](https://github.com/kokkos/kokkos/pull/3362)
- Kokkos Graph: add Cuda Graph implementation [\#3369](https://github.com/kokkos/kokkos/pull/3369)
- TeamPolicy: implemented autotuning of team sizes and vector lengths [\#3206](https://github.com/kokkos/kokkos/pull/3206)
- RangePolicy: Initialize all data members in default constructor [\#3509](https://github.com/kokkos/kokkos/pull/3509)

**Implemented enhancements BuildSystem:**
- Auto-generate core test files for all backends [\#3488](https://github.com/kokkos/kokkos/pull/3488)
- Avoid rewriting test files when calling cmake [\#3548](https://github.com/kokkos/kokkos/pull/3548)
- RULE\_LAUNCH\_COMPILE and RULE\_LAUNCH\_LINK system for nvcc\_wrapper [\#3136](https://github.com/kokkos/kokkos/pull/3136)
- Adding -include as a known argument to nvcc\_wrapper [\#3434](https://github.com/kokkos/kokkos/pull/3434)
- Install hpcbind script [\#3402](https://github.com/kokkos/kokkos/pull/3402)
- cmake/kokkos\_tribits.cmake: add parsing for args [\#3457](https://github.com/kokkos/kokkos/pull/3457)

**Implemented enhancements Tools:**
- Changed namespacing of Kokkos::Tools::Impl::Impl::tune\_policy [\#3455](https://github.com/kokkos/kokkos/pull/3455)
- Delegate to an impl allocate/deallocate method to allow specifying a SpaceHandle for MemorySpaces [\#3530](https://github.com/kokkos/kokkos/pull/3530)
- Use the Kokkos Profiling interface rather than the Impl interface [\#3518](https://github.com/kokkos/kokkos/pull/3518)
- Runtime option for tuning [\#3459](https://github.com/kokkos/kokkos/pull/3459)
- Dual View Tool Events [\#3326](https://github.com/kokkos/kokkos/pull/3326)

**Implemented enhancements Other:**
- Abort on errors instead of just printing [\#3528](https://github.com/kokkos/kokkos/pull/3528)
- Enable C++14 macros unconditionally [\#3449](https://github.com/kokkos/kokkos/pull/3449)
- Make ViewMapping trivially copyable [\#3436](https://github.com/kokkos/kokkos/pull/3436)
- Rename struct ViewMapping to class [\#3435](https://github.com/kokkos/kokkos/pull/3435)
- Replace enums in Kokkos\_ViewMapping.hpp (removes -Wextra) [\#3422](https://github.com/kokkos/kokkos/pull/3422)
- Use bool for enums representing bools [\#3416](https://github.com/kokkos/kokkos/pull/3416)
- Fence active instead of default execution space instances [\#3388](https://github.com/kokkos/kokkos/pull/3388)
- Refactor parallel\_reduce fence usage [\#3359](https://github.com/kokkos/kokkos/pull/3359)
- Moved Space EBO helpers to Kokkos\_EBO [\#3357](https://github.com/kokkos/kokkos/pull/3357)
- Add remove\_cvref type trait [\#3340](https://github.com/kokkos/kokkos/pull/3340)
- Adding identity type traits and update definition of identity\_t alias [\#3339](https://github.com/kokkos/kokkos/pull/3339)
- Add is\_specialization\_of type trait [\#3338](https://github.com/kokkos/kokkos/pull/3338)
- Make ScratchMemorySpace semi-regular [\#3309](https://github.com/kokkos/kokkos/pull/3309)
- Optimize min/max atomics with early exit on no-op case [\#3265](https://github.com/kokkos/kokkos/pull/3265)
- Refactor Backend Development [\#2941](https://github.com/kokkos/kokkos/pull/2941)

**Fixed bugs:**
- Fixup MDRangePolicy construction from Kokkos arrays [\#3591](https://github.com/kokkos/kokkos/pull/3591)
- Add atomic functions for unsigned long long using gcc built-in [\#3588](https://github.com/kokkos/kokkos/pull/3588)
- Fixup silent pointless comparison with zero in checked\_narrow\_cast (compiler workaround) [\#3566](https://github.com/kokkos/kokkos/pull/3566)
- Fixes for ROCm 3.9 [\#3565](https://github.com/kokkos/kokkos/pull/3565)
- Fix windows build issues which crept in for the CUDA build [\#3532](https://github.com/kokkos/kokkos/pull/3532)
- HIP Fix atomics of large data types and clean up lock arrays [\#3529](https://github.com/kokkos/kokkos/pull/3529)
- Pthreads fix exception resulting from 0 grain size [\#3510](https://github.com/kokkos/kokkos/pull/3510)
- Fixup do not require atomic operation to be default constructible [\#3503](https://github.com/kokkos/kokkos/pull/3503)
- Fix race condition in HIP backend [\#3467](https://github.com/kokkos/kokkos/pull/3467)
- Replace KOKKOS\_DEBUG with KOKKOS\_ENABLE\_DEBUG [\#3458](https://github.com/kokkos/kokkos/pull/3458)
- Fix multi-stream team scratch space definition for HIP [\#3398](https://github.com/kokkos/kokkos/pull/3398)
- HIP fix template deduction [\#3393](https://github.com/kokkos/kokkos/pull/3393)
- Fix compiling with HIP and C++17 [\#3390](https://github.com/kokkos/kokkos/pull/3390)
- Fix sigFPE in HIP blocksize deduction [\#3378](https://github.com/kokkos/kokkos/pull/3378)
- Type alias change: replace CS with CTS to avoid conflicts with NVSHMEM [\#3348](https://github.com/kokkos/kokkos/pull/3348)
- Clang compilation of CUDA backend on Windows [\#3345](https://github.com/kokkos/kokkos/pull/3345)
- Fix HBW support [\#3343](https://github.com/kokkos/kokkos/pull/3343)
- Added missing fences to unique token [\#3260](https://github.com/kokkos/kokkos/pull/3260)

**Incompatibilities:**
- Remove unused utilities (forward, move, and expand\_variadic) from Kokkos::Impl [\#3535](https://github.com/kokkos/kokkos/pull/3535)
- Remove unused traits [\#3534](https://github.com/kokkos/kokkos/pull/3534)
- HIP: Remove old HCC code [\#3301](https://github.com/kokkos/kokkos/pull/3301)
- Prepare for deprecation of ViewAllocateWithoutInitializing [\#3264](https://github.com/kokkos/kokkos/pull/3264)
- Remove ROCm backend [\#3148](https://github.com/kokkos/kokkos/pull/3148)

## [3.2.01](https://github.com/kokkos/kokkos/tree/3.2.01) (2020-11-17)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.2.00...3.2.01)

**Fixed bugs:**
- Disallow KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE in shared library builds [\#3332](https://github.com/kokkos/kokkos/pull/3332)
- Do not install libprinter-tool when testing is enabled [\#3313](https://github.com/kokkos/kokkos/pull/3313)
- Fix restrict/alignment following refactor [\#3373](https://github.com/kokkos/kokkos/pull/3373)
  - Intel fix: workaround compiler issue with using statement [\#3383](https://github.com/kokkos/kokkos/pull/3383)
- Fix zero-length reductions [#\3364](https://github.com/kokkos/kokkos/pull/3364)
  - Pthread zero-length reduction fix [\#3452](https://github.com/kokkos/kokkos/pull/3452)
  - HPX zero-length reduction fix [\#3470](https://github.com/kokkos/kokkos/pull/3470)
  - cuda/9.2 zero-length reduction fix [\#3580](https://github.com/kokkos/kokkos/pull/3580)
- Fix multi-stream scratch [#\3269](https://github.com/kokkos/kokkos/pull/3269)
- Guard KOKKOS_ALL_COMPILE_OPTIONS if Cuda is not enabled [\#3387](https://github.com/kokkos/kokkos/pull/3387)
- Do not include link flags for Fortran linkage [\#3384](https://github.com/kokkos/kokkos/pull/3384)
- Fix NVIDIA GPU arch macro with autodetection [\#3473](https://github.com/kokkos/kokkos/pull/3473)
- Fix libdl/test issues with Trilinos [\#3543](https://github.com/kokkos/kokkos/pull/3543)
  - Register Pthread as Tribits option to be enabled with Trilinos [\#3558](https://github.com/kokkos/kokkos/pull/3558)

**Implemented enhancements:**
- Separate Cuda timing-based tests into their own executable [\#3407](https://github.com/kokkos/kokkos/pull/3407)

## [3.2.00](https://github.com/kokkos/kokkos/tree/3.2.00) (2020-08-19)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.1.01...3.2.00)

**Implemented enhancements:**

- HIP:Enable stream in HIP [\#3163](https://github.com/kokkos/kokkos/issues/3163)
- HIP:Add support for shuffle reduction for the HIP backend [\#3154](https://github.com/kokkos/kokkos/issues/3154)
- HIP:Add implementations of missing HIPHostPinnedSpace methods for LAMMPS [\#3137](https://github.com/kokkos/kokkos/issues/3137)
- HIP:Require HIP 3.5.0 or higher [\#3099](https://github.com/kokkos/kokkos/issues/3099)
- HIP:WorkGraphPolicy for HIP [\#3096](https://github.com/kokkos/kokkos/issues/3096)
- OpenMPTarget: Significant update to the new experimental backend.  Requires C++17, works on Intel GPUs, reference counting fixes. [\#3169](https://github.com/kokkos/kokkos/issues/3169)
- Windows Cuda support [\#3018](https://github.com/kokkos/kokkos/issues/3018)
- Pass `-Wext-lambda-captures-this` to NVCC when support for `__host__ __device__` lambda is enabled from CUDA 11 [\#3241](https://github.com/kokkos/kokkos/issues/3241)
- Use explicit staging buffer for constant memory kernel launches and cleanup host/device synchronization [\#3234](https://github.com/kokkos/kokkos/issues/3234)
- Various fixup to policies including making TeamPolicy default constructible and making RangePolicy and TeamPolicy assignable: [\#3202](https://github.com/kokkos/kokkos/issues/3202) , [\#3203](https://github.com/kokkos/kokkos/issues/3203) , [\#3196](https://github.com/kokkos/kokkos/issues/3196)
- Annotations for `DefaultExectutionSpace` and `DefaultHostExectutionSpace` to use in static analysis [\#3189](https://github.com/kokkos/kokkos/issues/3189)
- Add documentation on using Spack to install Kokkos and developing packages that depend on Kokkos [\#3187](https://github.com/kokkos/kokkos/issues/3187)
- Add OpenMPTarget backend flags for NVC++ compiler [\#3185](https://github.com/kokkos/kokkos/issues/3185)
- Move deep\_copy/create\_mirror\_view on Experimental::OffsetView into Kokkos:: namespace [\#3166](https://github.com/kokkos/kokkos/issues/3166)
- Allow for larger block size in HIP [\#3165](https://github.com/kokkos/kokkos/issues/3165)
- View: Added names of Views to the different View initialize/free kernels [\#3159](https://github.com/kokkos/kokkos/issues/3159)
- Cuda: Caching cudaFunctorAttributes and whether L1/Shmem prefer was set [\#3151](https://github.com/kokkos/kokkos/issues/3151)
- BuildSystem: Improved performance in default configuration by defaulting to Release build [\#3131](https://github.com/kokkos/kokkos/issues/3131)
- Cuda: Update CUDA occupancy calculation [\#3124](https://github.com/kokkos/kokkos/issues/3124)
- Vector: Adding data() to Vector [\#3123](https://github.com/kokkos/kokkos/issues/3123)
- BuildSystem: Add CUDA Ampere configuration support [\#3122](https://github.com/kokkos/kokkos/issues/3122)
- General: Apply [[noreturn]] to Kokkos::abort when applicable [\#3106](https://github.com/kokkos/kokkos/issues/3106)
- TeamPolicy: Validate storage level argument passed to TeamPolicy::set\_scratch\_size() [\#3098](https://github.com/kokkos/kokkos/issues/3098)
- BuildSystem: Make kokkos\_has\_string() function in Makefile.kokkos case insensitive [\#3091](https://github.com/kokkos/kokkos/issues/3091)
- Modify KOKKOS\_FUNCTION macro for clang-tidy analysis [\#3087](https://github.com/kokkos/kokkos/issues/3087)
- Move allocation profiling to allocate/deallocate calls [\#3084](https://github.com/kokkos/kokkos/issues/3084)
- BuildSystem: FATAL\_ERROR when attempting in-source build [\#3082](https://github.com/kokkos/kokkos/issues/3082)
- Change enums in ScatterView to types [\#3076](https://github.com/kokkos/kokkos/issues/3076)
- HIP: Changes for new compiler/runtime [\#3067](https://github.com/kokkos/kokkos/issues/3067)
- Extract and use get\_gpu [\#3061](https://github.com/kokkos/kokkos/issues/3061) , [\#3048](https://github.com/kokkos/kokkos/issues/3048)
- Add is\_allocated to View-like containers [\#3059](https://github.com/kokkos/kokkos/issues/3059)
- Combined reducers for scalar references [\#3052](https://github.com/kokkos/kokkos/issues/3052)
- Add configurable capacity for UniqueToken [\#3051](https://github.com/kokkos/kokkos/issues/3051)
- Add installation testing [\#3034](https://github.com/kokkos/kokkos/issues/3034)
- HIP: Add UniqueToken [\#3020](https://github.com/kokkos/kokkos/issues/3020)
- Autodetect number of devices [\#3013](https://github.com/kokkos/kokkos/issues/3013)


**Fixed bugs:**

- Check error code from `cudaStreamSynchronize` in CUDA fences [\#3255](https://github.com/kokkos/kokkos/issues/3255)
- Fix issue with C++ standard flags when using `nvcc\_wrapper` with PGI [\#3254](https://github.com/kokkos/kokkos/issues/3254)
- Add missing threadfence in lock-based atomics [\#3208](https://github.com/kokkos/kokkos/issues/3208)
- Fix dedup of linker flags for shared lib on CMake <=3.12 [\#3176](https://github.com/kokkos/kokkos/issues/3176)
- Fix memory leak with CUDA streams [\#3170](https://github.com/kokkos/kokkos/issues/3170)
- BuildSystem: Fix OpenMP Target flags for Cray [\#3161](https://github.com/kokkos/kokkos/issues/3161)
- ScatterView: fix for OpenmpTarget remove inheritance from reducers [\#3162](https://github.com/kokkos/kokkos/issues/3162)
- BuildSystem: Set OpenMP flags according to host compiler [\#3127](https://github.com/kokkos/kokkos/issues/3127)
- OpenMP: Fix logic for nested omp in partition\_master bug [\#3101](https://github.com/kokkos/kokkos/issues/3101)
- nvcc\_wrapper: send --cudart to nvcc instead of host compiler [\#3092](https://github.com/kokkos/kokkos/issues/3092)
- BuildSystem: Fixes for Cuda/11 and c++17 [\#3085](https://github.com/kokkos/kokkos/issues/3085)
- HIP: Fix print\_configuration [\#3080](https://github.com/kokkos/kokkos/issues/3080)
- Conditionally define get\_gpu [\#3072](https://github.com/kokkos/kokkos/issues/3072)
- Fix bounds for ranges in random number generator [\#3069](https://github.com/kokkos/kokkos/issues/3069)
- Fix Cuda minor arch check [\#3035](https://github.com/kokkos/kokkos/issues/3035)
- BuildSystem: Add -expt-relaxed-constexpr flag to nvcc\_wrapper [\#3021](https://github.com/kokkos/kokkos/issues/3021)

**Incompatibilities:**

- Remove ETI support [\#3157](https://github.com/kokkos/kokkos/issues/3157)
- Remove KOKKOS\_INTERNAL\_ENABLE\_NON\_CUDA\_BACKEND [\#3147](https://github.com/kokkos/kokkos/issues/3147)
- Remove core/unit\_test/config [\#3146](https://github.com/kokkos/kokkos/issues/3146)
- Removed the preprocessor branch for KOKKOS\_ENABLE\_PROFILING [\#3115](https://github.com/kokkos/kokkos/issues/3115)
- Disable profiling with MSVC [\#3066](https://github.com/kokkos/kokkos/issues/3066)

**Closed issues:**

- Silent error (Validate storage level arg to set_scratch_size) [\#3097](https://github.com/kokkos/kokkos/issues/3097)
- Remove KOKKKOS\_ENABLE\_PROFILING Option [\#3095](https://github.com/kokkos/kokkos/issues/3095)
- Cuda 11 -\> allow C++17 [\#3083](https://github.com/kokkos/kokkos/issues/3083)
- In source build failure not explained [\#3081](https://github.com/kokkos/kokkos/issues/3081)
- Allow naming of Views for initialization kernel [\#3070](https://github.com/kokkos/kokkos/issues/3070)
- DefaultInit tests failing when using CTest resource allocation feature [\#3040](https://github.com/kokkos/kokkos/issues/3040)
- Add installation testing.  [\#3037](https://github.com/kokkos/kokkos/issues/3037)
- nvcc\_wrapper needs to handle `-expt-relaxed-constexpr` flag [\#3017](https://github.com/kokkos/kokkos/issues/3017)
- CPU core oversubscription warning on macOS with OpenMP backend [\#2996](https://github.com/kokkos/kokkos/issues/2996)
- Default behavior of KOKKOS\_NUM\_DEVICES to use all devices available [\#2975](https://github.com/kokkos/kokkos/issues/2975)
- Assert blocksize \> 0 [\#2974](https://github.com/kokkos/kokkos/issues/2974)
- Add ability to assign kokkos profile function from executable  [\#2973](https://github.com/kokkos/kokkos/issues/2973)
- ScatterView Support for the pre/post increment operator [\#2967](https://github.com/kokkos/kokkos/issues/2967)

- Compiler issue: Cuda build with clang 10 has errors with the atomic unit tests [\#3237](https://github.com/kokkos/kokkos/issues/3237)
- Incompatibility of flags for C++ standard with PGI v20.4 on Power9/NVIDIA V100 system [\#3252](https://github.com/kokkos/kokkos/issues/3252)
- Error configuring as subproject [\#3140](https://github.com/kokkos/kokkos/issues/3140)
- CMake fails with Nvidia compilers when the GPU architecture option is not supplied (Fix configure with OMPT and Cuda) [\#3207](https://github.com/kokkos/kokkos/issues/3207)
- PGI compiler being passed the gcc -fopenmp flag [\#3125](https://github.com/kokkos/kokkos/issues/3125)
- Cuda: Memory leak when using CUDA stream [\#3167](https://github.com/kokkos/kokkos/issues/3167)
- RangePolicy has an implicitly deleted assignment operator [\#3192](https://github.com/kokkos/kokkos/issues/3192)
- MemorySpace::allocate needs to have memory pool counting.  [\#3064](https://github.com/kokkos/kokkos/issues/3064)
- Missing write fence for lock based atomics on CUDA [\#3038](https://github.com/kokkos/kokkos/issues/3038)
- CUDA compute capability version check problem [\#3026](https://github.com/kokkos/kokkos/issues/3026)
- Make DynRankView fencing consistent [\#3014](https://github.com/kokkos/kokkos/issues/3014)
- nvcc\_wrapper cant handle -Xcompiler -o out.o [\#2993](https://github.com/kokkos/kokkos/issues/2993)
- Reductions of non-trivial types of size 4 fail in CUDA shfl operations [\#2990](https://github.com/kokkos/kokkos/issues/2990)
- complex\_double misalignment in reduce, clang+CUDA [\#2989](https://github.com/kokkos/kokkos/issues/2989)
- Span of degenerated \(zero-length\) subviews is not zero in some special cases [\#2979](https://github.com/kokkos/kokkos/issues/2979)
- Rank 1 custom layouts dont work as expected. [\#2840](https://github.com/kokkos/kokkos/issues/2840)

## [3.1.01](https://github.com/kokkos/kokkos/tree/3.1.1) (2020-04-14)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.1.00...3.1.1)

**Fixed bugs:**

- Fix complex_double misalignment in reduce, clang+CUDA [\#2989](https://github.com/kokkos/kokkos/issues/2989)
- Fix compilation fails when profiling disabled and CUDA enabled [\#3001](https://github.com/kokkos/kokkos/issues/3001)
- Fix cuda reduction of non-trivial scalars of size 4 [\#2990](https://github.com/kokkos/kokkos/issues/2990)
- Configure and install version file when building in Trilinos [\#2957](https://github.com/kokkos/kokkos/pull/2957)
- Fix OpenMPTarget build missing include and namespace [\#3000](https://github.com/kokkos/kokkos/issues/3000)
- fix typo in KOKKOS_SET_EXE_PROPERTY() [\#2959](https://github.com/kokkos/kokkos/issues/2959)
- Fix non-zero span subviews of zero sized subviews [\#2979](https://github.com/kokkos/kokkos/issues/2979)

## [3.1.00](https://github.com/kokkos/kokkos/tree/3.1.00) (2020-04-14)
[Full Changelog](https://github.com/kokkos/kokkos/compare/3.0.00...3.1.00)

**Features:**

- HIP Support for AMD
- OpenMPTarget Support with clang
- Windows VS19 (Serial) Support [\#1533](https://github.com/kokkos/kokkos/issues/1533)

**Implemented enhancements:**

- generate\_makefile.bash should allow tests to be disabled [\#2886](https://github.com/kokkos/kokkos/issues/2886)
- clang/7+cuda/9 build -Werror-unused parameter error in nightly test [\#2884](https://github.com/kokkos/kokkos/issues/2884)
- ScatterView memory space is not user settable [\#2826](https://github.com/kokkos/kokkos/issues/2826)
- clang/8+cuda/10.0 build error with c++17 [\#2809](https://github.com/kokkos/kokkos/issues/2809)
- warnings.... [\#2805](https://github.com/kokkos/kokkos/issues/2805)
- Kokkos version in cpp define [\#2787](https://github.com/kokkos/kokkos/issues/2787)
- Remove Defunct QThreads Backend [\#2751](https://github.com/kokkos/kokkos/issues/2751)
- Improve Kokkos::fence behavior with multiple execution spaces [\#2659](https://github.com/kokkos/kokkos/issues/2659)
- polylithic\(?\) initialization of Kokkos [\#2658](https://github.com/kokkos/kokkos/issues/2658)
- Unnecessary\(?\) check for host execution space initialization from Cuda initialization [\#2652](https://github.com/kokkos/kokkos/issues/2652)
- Kokkos error reporting failures with CUDA GPUs in exclusive mode [\#2471](https://github.com/kokkos/kokkos/issues/2471)
- atomicMax equivalent \(and other atomics\) [\#2401](https://github.com/kokkos/kokkos/issues/2401)
- Fix alignment for Kokkos::complex [\#2255](https://github.com/kokkos/kokkos/issues/2255)
- Warnings with Cuda 10.1 [\#2206](https://github.com/kokkos/kokkos/issues/2206)
- dual view with Kokkos::ViewAllocateWithoutInitializing [\#2188](https://github.com/kokkos/kokkos/issues/2188)
- Check error code  from cudaOccupancyMaxActiveBlocksPerMultiprocessor [\#2172](https://github.com/kokkos/kokkos/issues/2172)
- Add non-member Kokkos::resize/realloc for DualView [\#2170](https://github.com/kokkos/kokkos/issues/2170)
- Construct DualView without initialization [\#2046](https://github.com/kokkos/kokkos/issues/2046)
- Expose is\_assignable to determine if one view can be assigned to another [\#1936](https://github.com/kokkos/kokkos/issues/1936)
- profiling label [\#1935](https://github.com/kokkos/kokkos/issues/1935)
- team\_broadcast of bool failed on CUDA backend [\#1908](https://github.com/kokkos/kokkos/issues/1908)
- View static\_extent [\#660](https://github.com/kokkos/kokkos/issues/660)
- Misleading Kokkos::Cuda::initialize ERROR message when compiled for wrong GPU architecture [\#1944](https://github.com/kokkos/kokkos/issues/1944)
- Cryptic Error When Malloc Fails [\#2164](https://github.com/kokkos/kokkos/issues/2164)
- Drop support for intermediate standards in CMake [\#2336](https://github.com/kokkos/kokkos/issues/2336)

**Fixed bugs:**

- DualView sync\_device with length zero creates cuda errors [\#2946](https://github.com/kokkos/kokkos/issues/2946)
- building with nvcc and clang \(or clang based XL\) as host compiler: "Kokkos::atomic\_fetch\_min\(volatile int \*, int\)" has already been defined [\#2903](https://github.com/kokkos/kokkos/issues/2903)
- Cuda 9.1,10.1 debug builds failing due to -Werror=unused-parameter [\#2880](https://github.com/kokkos/kokkos/issues/2880)
- clang -Werror: Kokkos\_FixedBufferMemoryPool.hpp:140:28: error: unused parameter 'alloc\_size' [\#2869](https://github.com/kokkos/kokkos/issues/2869)
- intel/16.0.1, intel/17.0.1 nightly build failures with debugging enabled [\#2867](https://github.com/kokkos/kokkos/issues/2867)
- intel/16.0.1 debug build errors [\#2863](https://github.com/kokkos/kokkos/issues/2863)
- xl/16.1.1 with cpp14, openmp build, nightly test failures [\#2856](https://github.com/kokkos/kokkos/issues/2856)
- Intel nightly test failures: team\_vector [\#2852](https://github.com/kokkos/kokkos/issues/2852)
- Kokkos Views with intmax/2\<N\<intmax can hang during construction [\#2850](https://github.com/kokkos/kokkos/issues/2850)
- workgraph\_fib test seg-faults with threads backend and hwloc [\#2797](https://github.com/kokkos/kokkos/issues/2797)
- cuda.view\_64bit test hangs on Power8+Kepler37 system - develop and 2.9.00 branches [\#2771](https://github.com/kokkos/kokkos/issues/2771)
- device\_type for Kokkos\_Random ?  [\#2693](https://github.com/kokkos/kokkos/issues/2693)
- "More than one tag given" error in Experimental::require\(\) [\#2608](https://github.com/kokkos/kokkos/issues/2608)
- Segfault on Marvell from our finalization stack [\#2542](https://github.com/kokkos/kokkos/issues/2542)

## [3.0.00](https://github.com/kokkos/kokkos/tree/3.0.00) (2020-01-27)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.9.00...3.0.00)

**Implemented enhancements:**

- BuildSystem: Standalone Modern CMake Support [\#2104](https://github.com/kokkos/kokkos/issues/2104)
- StyleFormat: ClangFormat Style [\#2157](https://github.com/kokkos/kokkos/issues/2157)
- Documentation: Document build system and CMake philosophy [\#2263](https://github.com/kokkos/kokkos/issues/2263)
- BuildSystem: Add Alias with Namespace Kokkos:: to Interal Libraries [\#2530](https://github.com/kokkos/kokkos/issues/2530)
- BuildSystem: Universal Kokkos find\_package [\#2099](https://github.com/kokkos/kokkos/issues/2099)
- BuildSystem: Dropping support for Kokkos\_{DEVICES,OPTIONS,ARCH} in CMake [\#2329](https://github.com/kokkos/kokkos/issues/2329)
- BuildSystem: Set Kokkos\_DEVICES and Kokkos\_ARCH variables in exported CMake configuration [\#2193](https://github.com/kokkos/kokkos/issues/2193)
- BuildSystem: Drop support for CUDA 7 and CUDA 8 [\#2489](https://github.com/kokkos/kokkos/issues/2489)
- BuildSystem: Drop CMake option SEPARATE\_TESTS [\#2266](https://github.com/kokkos/kokkos/issues/2266)
- BuildSystem: Support expt-relaxed-constexpr same as expt-extended-lambda [\#2411](https://github.com/kokkos/kokkos/issues/2411)
- BuildSystem: Add Xnvlink to command line options allowed in nvcc\_wrapper [\#2197](https://github.com/kokkos/kokkos/issues/2197)
- BuildSystem: Install Kokkos config files and target files to lib/cmake/Kokkos [\#2162](https://github.com/kokkos/kokkos/issues/2162)
- BuildSystem: nvcc\_wrappers and c++ 14 [\#2035](https://github.com/kokkos/kokkos/issues/2035)
- BuildSystem: Kokkos version major/version minor \(Feature request\) [\#1930](https://github.com/kokkos/kokkos/issues/1930)
- BuildSystem: CMake namespaces \(and other modern cmake cleanup\) [\#1924](https://github.com/kokkos/kokkos/issues/1924)
- BuildSystem: Remove capability to install Kokkos via GNU Makefiles [\#2332](https://github.com/kokkos/kokkos/issues/2332)
- Documentation: Remove PDF ProgrammingGuide in Kokkos replace with link [\#2244](https://github.com/kokkos/kokkos/issues/2244)
- View: Add Method to Resize View without Initialization [\#2048](https://github.com/kokkos/kokkos/issues/2048)
- Vector: implement “insert” method for Kokkos\_Vector  \(as a serial function on host\) [\#2437](https://github.com/kokkos/kokkos/issues/2437)

**Fixed bugs:**

- ParallelScan: Kokkos::parallel\scan fix race condition seen in inter-block fence [\#2681](https://github.com/kokkos/kokkos/issues/2681)
- OffsetView: Kokkos::OffsetView missing constructor which takes pointer [\#2247](https://github.com/kokkos/kokkos/issues/2247)
- OffsetView: Kokkos::OffsetView: allow offset=0 [\#2246](https://github.com/kokkos/kokkos/issues/2246)
- DeepCopy: Missing DeepCopy instrumentation in Kokkos [\#2522](https://github.com/kokkos/kokkos/issues/2522)
- nvcc\_wrapper: --host-only fails with multiple -W\* flags [\#2484](https://github.com/kokkos/kokkos/issues/2484)
- nvcc\_wrapper: taking first -std option is counterintuitive [\#2553](https://github.com/kokkos/kokkos/issues/2553)
- Subview: Error taking subviews of views with static_extents of min rank [\#2448](https://github.com/kokkos/kokkos/issues/2448)
- TeamPolicy: reducers with valuetypes without += broken on CUDA [\#2410](https://github.com/kokkos/kokkos/issues/2410)
- Libs: Fix inconsistency of Kokkos library names in Kokkos and Trilinos [\#1902](https://github.com/kokkos/kokkos/issues/1902)
- Complex: operator\>\> for complex\<T\> uses std::ostream, not std::istream [\#2313](https://github.com/kokkos/kokkos/issues/2313)
- Macros: Restrict not honored for non-intel compilers  [\#1922](https://github.com/kokkos/kokkos/issues/1922)


## [2.9.00](https://github.com/kokkos/kokkos/tree/2.9.00) (2019-06-24)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.8.00...2.9.00)

**Implemented enhancements:**

- Capability: CUDA Streams [\#1723](https://github.com/kokkos/kokkos/issues/1723)
- Capability: CUDA Stream support for parallel\_reduce [\#2061](https://github.com/kokkos/kokkos/issues/2061)
- Capability: Feature Request: TeamVectorRange [\#713](https://github.com/kokkos/kokkos/issues/713)
- Capability: Adding HPX backend [\#2080](https://github.com/kokkos/kokkos/issues/2080)
- Capability: TaskScheduler to have multiple queues [\#565](https://github.com/kokkos/kokkos/issues/565)
- Capability: Support for additional reductions in ScatterView [\#1674](https://github.com/kokkos/kokkos/issues/1674)
- Capability: Request: deep\_copy within parallel regions [\#689](https://github.com/kokkos/kokkos/issues/689)
- Capability: Feature Request: `create\_mirror\_view\_without\_initializing` [\#1765](https://github.com/kokkos/kokkos/issues/1765)
- View: Use SFINAE to restrict possible View type conversions [\#2127](https://github.com/kokkos/kokkos/issues/2127)
- Deprecation: Deprecate ExecutionSpace::fence\(\) as static function and make it non-static [\#2140](https://github.com/kokkos/kokkos/issues/2140)
- Deprecation: Deprecate LayoutTileLeft [\#2122](https://github.com/kokkos/kokkos/issues/2122)
- Macros: KOKKOS\_RESTRICT defined for non-Intel compilers [\#2038](https://github.com/kokkos/kokkos/issues/2038)

**Fixed bugs:**

- Cuda: TeamThreadRange loop count on device is passed by reference to host static constexpr [\#1733](https://github.com/kokkos/kokkos/issues/1733)
- Cuda: Build error with relocatable device code with CUDA 10.1 GCC 7.3 [\#2134](https://github.com/kokkos/kokkos/issues/2134)
- Cuda: cudaFuncSetCacheConfig is setting CachePreferShared too often [\#2066](https://github.com/kokkos/kokkos/issues/2066)
- Cuda: TeamPolicy doesn't throw then created with non-viable vector length and also doesn't backscale to viable one [\#2020](https://github.com/kokkos/kokkos/issues/2020)
- Cuda: cudaMemcpy error for large league sizes on V100 [\#1991](https://github.com/kokkos/kokkos/issues/1991)
- Cuda: illegal warp sync in parallel\_reduce by functor on Turing 75 [\#1958](https://github.com/kokkos/kokkos/issues/1958)
- TeamThreadRange: Inconsistent results from TeamThreadRange reduction [\#1905](https://github.com/kokkos/kokkos/issues/1905)
- Atomics: atomic\_fetch\_oper & atomic\_oper\_fetch don't build for complex\<float\> [\#1964](https://github.com/kokkos/kokkos/issues/1964)
- Views: Kokkos randomread Views leak memory [\#2155](https://github.com/kokkos/kokkos/issues/2155)
- ScatterView: LayoutLeft overload currently non-functional [\#2165](https://github.com/kokkos/kokkos/issues/2165)
- KNL: With intel 17.2.174 illegal instruction in random number test [\#2078](https://github.com/kokkos/kokkos/issues/2078)
- Bitset: Enable copy constructor on device [\#2094](https://github.com/kokkos/kokkos/issues/2094)
- Examples: do not compile due to template deduction error \(multi\_fem\) [\#1928](https://github.com/kokkos/kokkos/issues/1928)

## [2.8.00](https://github.com/kokkos/kokkos/tree/2.8.00) (2019-02-05)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.7.24...2.8.00)

**Implemented enhancements:**

- Capability, Tests: C++14 support and testing [\#1914](https://github.com/kokkos/kokkos/issues/1914)
- Capability: Add environment variables for all command line arguments [\#1798](https://github.com/kokkos/kokkos/issues/1798)
- Capability: --kokkos-ndevices not working for Slurm [\#1920](https://github.com/kokkos/kokkos/issues/1920)
- View: Undefined behavior when deep copying from and to an empty unmanaged view [\#1967](https://github.com/kokkos/kokkos/issues/1967)
- BuildSystem: nvcc\_wrapper should stop immediately if nvcc is not in PATH [\#1861](https://github.com/kokkos/kokkos/issues/1861)

**Fixed bugs:**

- Cuda: Fix Volta Issues 1 Non-deterministic behavior on Volta, runs fine on Pascal [\#1949](https://github.com/kokkos/kokkos/issues/1949)
- Cuda: Fix Volta Issues 2 CUDA Team Scan gives wrong values on Volta with -G compile flag [\#1942](https://github.com/kokkos/kokkos/issues/1942)
- Cuda: illegal warp sync in parallel\_reduce by functor on Turing 75 [\#1958](https://github.com/kokkos/kokkos/issues/1958)
- Threads: Pthreads backend does not handle RangePolicy with offset correctly [\#1976](https://github.com/kokkos/kokkos/issues/1976)
- Atomics: atomic\_fetch\_oper has no case for Kokkos::complex\<double\> or other 16-byte types [\#1951](https://github.com/kokkos/kokkos/issues/1951)
- MDRangePolicy: Fix zero-length range [\#1948](https://github.com/kokkos/kokkos/issues/1948)
- TeamThreadRange: TeamThreadRange MaxLoc reduce doesnt compile  [\#1909](https://github.com/kokkos/kokkos/issues/1909)

## [2.7.24](https://github.com/kokkos/kokkos/tree/2.7.24) (2018-11-04)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.7.00...2.7.24)

**Implemented enhancements:**

- DualView: Add non-templated functions for sync, need\_sync, view, modify [\#1858](https://github.com/kokkos/kokkos/issues/1858)
- DualView: Avoid needlessly allocates and initializes modify\_host and modify\_device flag views [\#1831](https://github.com/kokkos/kokkos/issues/1831)
- DualView: Incorrect deduction of "not device type" [\#1659](https://github.com/kokkos/kokkos/issues/1659)
- BuildSystem: Add KOKKOS\_ENABLE\_CXX14 and KOKKOS\_ENABLE\_CXX17 [\#1602](https://github.com/kokkos/kokkos/issues/1602)
- BuildSystem: Installed kokkos\_generated\_settings.cmake contains build directories instead of install directories [\#1838](https://github.com/kokkos/kokkos/issues/1838)
- BuildSystem: KOKKOS\_ARCH: add ticks to printout of improper arch setting [\#1649](https://github.com/kokkos/kokkos/issues/1649)
- BuildSystem: Make core/src/Makefile for Cuda use needed nvcc\_wrapper [\#1296](https://github.com/kokkos/kokkos/issues/1296)
- Build: Support PGI as host compiler for NVCC [\#1828](https://github.com/kokkos/kokkos/issues/1828)
- Build: Many Warnings Fixed e.g.[\#1786](https://github.com/kokkos/kokkos/issues/1786)
- Capability: OffsetView with non-zero begin index [\#567](https://github.com/kokkos/kokkos/issues/567)
- Capability: Reductions into device side view [\#1788](https://github.com/kokkos/kokkos/issues/1788)
- Capability: Add max\_size to Kokkos::Array [\#1760](https://github.com/kokkos/kokkos/issues/1760)
- Capability: View Assignment: LayoutStride -\> LayoutLeft and LayoutStride -\> LayoutRight [\#1594](https://github.com/kokkos/kokkos/issues/1594)
- Capability: Atomic function allow implicit conversion of update argument [\#1571](https://github.com/kokkos/kokkos/issues/1571)
- Capability: Add team\_size\_max with tagged functors [\#663](https://github.com/kokkos/kokkos/issues/663)
- Capability: Fix allignment of views from Kokkos\_ScratchSpace should use different alignment [\#1700](https://github.com/kokkos/kokkos/issues/1700)
- Capabilitiy: create\_mirror\_view\_and\_copy for DynRankView [\#1651](https://github.com/kokkos/kokkos/issues/1651)
- Capability: DeepCopy HBWSpace / HostSpace [\#548](https://github.com/kokkos/kokkos/issues/548)
- ROCm: support team vector scan  [\#1645](https://github.com/kokkos/kokkos/issues/1645)
- ROCm:  Merge from rocm-hackathon2 [\#1636](https://github.com/kokkos/kokkos/issues/1636)
- ROCm:  Add ParallelScanWithTotal [\#1611](https://github.com/kokkos/kokkos/issues/1611)
- ROCm: Implement MDRange in ROCm [\#1314](https://github.com/kokkos/kokkos/issues/1314)
- ROCm: Implement Reducers for Nested Parallelism Levels [\#963](https://github.com/kokkos/kokkos/issues/963)
- ROCm: Add asynchronous deep copy [\#959](https://github.com/kokkos/kokkos/issues/959)
- Tests: Memory pool test seems to allocate 8GB [\#1830](https://github.com/kokkos/kokkos/issues/1830)
- Tests: Add unit\_test for team\_broadcast [\#734](https://github.com/kokkos/kokkos/issues/734)

**Fixed bugs:**

- BuildSystem: Makefile.kokkos gets gcc-toolchain wrong if gcc is cached [\#1841](https://github.com/kokkos/kokkos/issues/1841)
- BuildSystem: kokkos\_generated\_settings.cmake placement is inconsistent [\#1771](https://github.com/kokkos/kokkos/issues/1771)
- BuildSystem: Invalid escape sequence \. in kokkos\_functions.cmake [\#1661](https://github.com/kokkos/kokkos/issues/1661)
- BuildSystem: Problem in Kokkos generated cmake file [\#1770](https://github.com/kokkos/kokkos/issues/1770)
- BuildSystem: invalid file names on windows [\#1671](https://github.com/kokkos/kokkos/issues/1671)
- Tests: reducers min/max\_loc test fails randomly due to multiple min values and thus multiple valid locations [\#1681](https://github.com/kokkos/kokkos/issues/1681)
- Tests: cuda.scatterview unit test causes "Bus error" when force\_uvm and enable\_lambda are enabled [\#1852](https://github.com/kokkos/kokkos/issues/1852)
- Tests: cuda.cxx11 unit test fails when force\_uvm and enable\_lambda are enabled [\#1850](https://github.com/kokkos/kokkos/issues/1850)
- Tests: threads.reduce\_device\_view\_range\_policy failing with Cuda/8.0.44 and RDC [\#1836](https://github.com/kokkos/kokkos/issues/1836)
- Build: compile error when compiling Kokkos with hwloc 2.0.1 \(on OSX 10.12.6, with g++ 7.2.0\) [\#1506](https://github.com/kokkos/kokkos/issues/1506)
- Build: dual\_view.view broken with UVM [\#1834](https://github.com/kokkos/kokkos/issues/1834)
- Build: White cuda/9.2 + gcc/7.2 warnings triggering errors  [\#1833](https://github.com/kokkos/kokkos/issues/1833)
- Build: warning: enum constant in boolean context [\#1813](https://github.com/kokkos/kokkos/issues/1813)
- Capability: Fix overly conservative max\_team\_size thingy [\#1808](https://github.com/kokkos/kokkos/issues/1808)
- DynRankView: Ctors taking ViewAllocateWithoutInitializing broken [\#1783](https://github.com/kokkos/kokkos/issues/1783)
- Cuda: Apollo cuda.team\_broadcast test fail with clang-6.0 [\#1762](https://github.com/kokkos/kokkos/issues/1762)
- Cuda: Clang spurious test failure in impl\_view\_accessible [\#1753](https://github.com/kokkos/kokkos/issues/1753)
- Cuda: Kokkos::complex\<double\> atomic deadlocks with Clang 6 Cuda build with -O0 [\#1752](https://github.com/kokkos/kokkos/issues/1752)
- Cuda: LayoutStride Test fails for UVM as default memory space [\#1688](https://github.com/kokkos/kokkos/issues/1688)
- Cuda: Scan wrong values on Volta [\#1676](https://github.com/kokkos/kokkos/issues/1676)
- Cuda: Kokkos::deep\_copy error with CudaUVM and Kokkos::Serial spaces [\#1652](https://github.com/kokkos/kokkos/issues/1652)
- Cuda: cudaErrorInvalidConfiguration with debug build [\#1647](https://github.com/kokkos/kokkos/issues/1647)
- Cuda: parallel\_for with TeamPolicy::team\_size\_recommended with launch bounds not working -- reported by Daniel Holladay [\#1283](https://github.com/kokkos/kokkos/issues/1283)
- Cuda: Using KOKKOS\_CLASS\_LAMBDA in a class with Kokkos::Random\_XorShift64\_Pool member data [\#1696](https://github.com/kokkos/kokkos/issues/1696)
- Long Build Times on Darwin [\#1721](https://github.com/kokkos/kokkos/issues/1721)
- Capability: Typo in Kokkos\_Sort.hpp - BinOp3D - wrong comparison [\#1720](https://github.com/kokkos/kokkos/issues/1720)
- Buffer overflow in SharedAllocationRecord in Kokkos\_HostSpace.cpp [\#1673](https://github.com/kokkos/kokkos/issues/1673)
- Serial unit test failure [\#1632](https://github.com/kokkos/kokkos/issues/1632)

## [2.7.00](https://github.com/kokkos/kokkos/tree/2.7.00) (2018-05-24)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.6.00...2.7.00)

**Part of the Kokkos C++ Performance Portability Programming EcoSystem 2.7**

**Implemented enhancements:**

- Deprecate team\_size auto adjusting to maximal value possible [\#1618](https://github.com/kokkos/kokkos/issues/1618)
- DynamicView - remove restrictions to std::is\_trivial types and value\_type is power of two [\#1586](https://github.com/kokkos/kokkos/issues/1586)
- Kokkos::StaticCrsGraph does not propagate memory traits \(e.g., Unmanaged\) [\#1581](https://github.com/kokkos/kokkos/issues/1581)
- Adding ETI for DeepCopy / ViewFill etc. [\#1578](https://github.com/kokkos/kokkos/issues/1578)
- Deprecate all the left over KOKKOS\_HAVE\_ Macros and Kokkos\_OldMacros.hpp [\#1572](https://github.com/kokkos/kokkos/issues/1572)
- Error if Kokkos\_ARCH set in CMake [\#1555](https://github.com/kokkos/kokkos/issues/1555)
- Deprecate ExecSpace::initialize / ExecSpace::finalize [\#1532](https://github.com/kokkos/kokkos/issues/1532)
- New API for TeamPolicy property setting [\#1531](https://github.com/kokkos/kokkos/issues/1531)
- clang 6.0 + cuda debug out-of-memory test failure [\#1521](https://github.com/kokkos/kokkos/issues/1521)
- Cuda UniqueToken interface not consistent with other backends [\#1505](https://github.com/kokkos/kokkos/issues/1505)
- Move Reducers out of Experimental namespace [\#1494](https://github.com/kokkos/kokkos/issues/1494)
- Provide scope guard for initialize/finalize [\#1479](https://github.com/kokkos/kokkos/issues/1479)
- Check Kokkos::is\_initialized in SharedAllocationRecord dtor [\#1465](https://github.com/kokkos/kokkos/issues/1465)
- Remove static list of allocations [\#1464](https://github.com/kokkos/kokkos/issues/1464)
- Makefiles: Support single compile/link line use case [\#1402](https://github.com/kokkos/kokkos/issues/1402)
- ThreadVectorRange with a range  [\#1400](https://github.com/kokkos/kokkos/issues/1400)
- Exclusive scan + last value API [\#1358](https://github.com/kokkos/kokkos/issues/1358)
- Install kokkos\_generated\_settings.cmake [\#1348](https://github.com/kokkos/kokkos/issues/1348)
- Kokkos arrays \(not views!\) don't do bounds checking in debug mode [\#1342](https://github.com/kokkos/kokkos/issues/1342)
- Expose round-robin GPU assignment outside of initialize\(int, char\*\*\) [\#1318](https://github.com/kokkos/kokkos/issues/1318)
- DynamicView misses use\_count and label function [\#1298](https://github.com/kokkos/kokkos/issues/1298)
- View constructor should check arguments [\#1286](https://github.com/kokkos/kokkos/issues/1286)
- False Positive on Oversubscription Warning [\#1207](https://github.com/kokkos/kokkos/issues/1207)
- Allow \(require\) execution space for 1st arg of VerifyExecutionCanAccessMemorySpace [\#1192](https://github.com/kokkos/kokkos/issues/1192)
- ROCm: Add ROCmHostPinnedSpace [\#958](https://github.com/kokkos/kokkos/issues/958)
- power of two functions [\#656](https://github.com/kokkos/kokkos/issues/656)
- CUDA 8 has 64bit \_\_shfl [\#361](https://github.com/kokkos/kokkos/issues/361)
- Add TriBITS/CMake configure information about node types [\#243](https://github.com/kokkos/kokkos/issues/243)

**Fixed bugs:**

- CUDA atomic\_fetch\_sub for doubles is hitting CAS instead of intrinsic [\#1624](https://github.com/kokkos/kokkos/issues/1624)
- Bug: use of ballot on Volta [\#1612](https://github.com/kokkos/kokkos/issues/1612)
- Kokkos::deep\_copy memory access failures [\#1583](https://github.com/kokkos/kokkos/issues/1583)
- g++ -std option doubly set for cmake project [\#1548](https://github.com/kokkos/kokkos/issues/1548)
- ViewFill for 1D Views of larger 32bit entries fails [\#1541](https://github.com/kokkos/kokkos/issues/1541)
- CUDA Volta another warpsync bug [\#1520](https://github.com/kokkos/kokkos/issues/1520)
- triple\_nested\_parallelism fails with KOKKOS\_DEBUG and CUDA [\#1513](https://github.com/kokkos/kokkos/issues/1513)
- Jenkins errors in Kokkos\_SharedAlloc.cpp with debug build [\#1511](https://github.com/kokkos/kokkos/issues/1511)
- Kokkos::Sort out-of-bounds with empty bins [\#1504](https://github.com/kokkos/kokkos/issues/1504)
- Get rid of deprecated functions inside Kokkos [\#1484](https://github.com/kokkos/kokkos/issues/1484)
- get\_work\_partition casts int64\_t to int, causing a seg fault [\#1481](https://github.com/kokkos/kokkos/issues/1481)
- NVCC bug with \_\_device\_\_ on defaulted function [\#1470](https://github.com/kokkos/kokkos/issues/1470)
- CMake example broken with CUDA backend [\#1468](https://github.com/kokkos/kokkos/issues/1468)


## [2.6.00](https://github.com/kokkos/kokkos/tree/2.6.00) (2018-03-07)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.5.00...2.6.00)

**Part of the Kokkos C++ Performance Portability Programming EcoSystem 2.6**

**Implemented enhancements:**

- Support NVIDIA Volta microarchitecture [\#1466](https://github.com/kokkos/kokkos/issues/1466)
- Kokkos - Define empty functions when profiling disabled [\#1424](https://github.com/kokkos/kokkos/issues/1424)
- Don't use \_\_constant\_\_ cache for lock arrays, enable once per run update instead of once per call [\#1385](https://github.com/kokkos/kokkos/issues/1385)
- task dag enhancement. [\#1354](https://github.com/kokkos/kokkos/issues/1354)
- Cuda task team collectives and stack size [\#1353](https://github.com/kokkos/kokkos/issues/1353)
- Replace View operator acceptance of more than rank integers with 'access' function [\#1333](https://github.com/kokkos/kokkos/issues/1333)
- Interoperability: Do not shut down backend execution space runtimes upon calling finalize. [\#1305](https://github.com/kokkos/kokkos/issues/1305)
- shmem\_size for LayoutStride [\#1291](https://github.com/kokkos/kokkos/issues/1291)
- Kokkos::resize performs poorly on 1D Views [\#1270](https://github.com/kokkos/kokkos/issues/1270)
- stride\(\) is inconsistent with dimension\(\), extent\(\), etc. [\#1214](https://github.com/kokkos/kokkos/issues/1214)
- Kokkos::sort defaults to std::sort on host [\#1208](https://github.com/kokkos/kokkos/issues/1208)
- DynamicView with host size grow [\#1206](https://github.com/kokkos/kokkos/issues/1206)
- Unmanaged View with Anonymous Memory Space [\#1175](https://github.com/kokkos/kokkos/issues/1175)
- Sort subset of Kokkos::DynamicView [\#1160](https://github.com/kokkos/kokkos/issues/1160)
- MDRange policy doesn't support lambda reductions [\#1054](https://github.com/kokkos/kokkos/issues/1054)
- Add ability to set hook on Kokkos::finalize [\#714](https://github.com/kokkos/kokkos/issues/714)
- Atomics with Serial Backend - Default should be Disable? [\#549](https://github.com/kokkos/kokkos/issues/549)
- KOKKOS\_ENABLE\_DEPRECATED\_CODE [\#1359](https://github.com/kokkos/kokkos/issues/1359)

**Fixed bugs:**

- cuda\_internal\_maximum\_warp\_count returns 8, but I believe it should return 16 for P100  [\#1269](https://github.com/kokkos/kokkos/issues/1269)
- Cuda: level 1 scratch memory bug \(reported by Stan Moore\) [\#1434](https://github.com/kokkos/kokkos/issues/1434)
- MDRangePolicy Reduction requires value\_type typedef in Functor [\#1379](https://github.com/kokkos/kokkos/issues/1379)
- Kokkos DeepCopy between empty views fails [\#1369](https://github.com/kokkos/kokkos/issues/1369)
- Several issues with new CMake build infrastructure \(reported by Eric Phipps\) [\#1365](https://github.com/kokkos/kokkos/issues/1365)
- deep\_copy between rank-1 host/device views of differing layouts without UVM no longer works \(reported by Eric Phipps\) [\#1363](https://github.com/kokkos/kokkos/issues/1363)
- Profiling can't be disabled in CMake, and a parallel\_for is missing for tasks \(reported by Kyungjoo Kim\) [\#1349](https://github.com/kokkos/kokkos/issues/1349)
- get\_work\_partition int overflow \(reported by berryj5\) [\#1327](https://github.com/kokkos/kokkos/issues/1327)
- Kokkos::deep\_copy must fence even if the two views are the same [\#1303](https://github.com/kokkos/kokkos/issues/1303)
- CudaUVMSpace::allocate/deallocate must fence [\#1302](https://github.com/kokkos/kokkos/issues/1302)
- ViewResize on CUDA fails in Debug because of too many resources requested [\#1299](https://github.com/kokkos/kokkos/issues/1299)
- Cuda 9 and intrepid2 calls from Panzer. [\#1183](https://github.com/kokkos/kokkos/issues/1183)
- Slowdown due to tracking\_enabled\(\) in 2.04.00 \(found by Albany app\) [\#1016](https://github.com/kokkos/kokkos/issues/1016)
- Bounds checking fails with zero-span Views \(reported by Stan Moore\) [\#1411](https://github.com/kokkos/kokkos/issues/1411)


## [2.5.00](https://github.com/kokkos/kokkos/tree/2.5.00) (2017-12-15)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.04.11...2.5.00)

**Part of the Kokkos C++ Performance Portability Programming EcoSystem 2.5**

**Implemented enhancements:**

- Provide Makefile.kokkos logic for CMake and TriBITS [\#878](https://github.com/kokkos/kokkos/issues/878)
- Add Scatter View [\#825](https://github.com/kokkos/kokkos/issues/825)
- Drop gcc 4.7 and intel 14 from supported compiler list [\#603](https://github.com/kokkos/kokkos/issues/603)
- Enable construction of unmanaged view using common\_view\_alloc\_prop [\#1170](https://github.com/kokkos/kokkos/issues/1170)
- Unused Function Warning with XL [\#1267](https://github.com/kokkos/kokkos/issues/1267)
- Add memory pool parameter check [\#1218](https://github.com/kokkos/kokkos/issues/1218)
- CUDA9: Fix warning for unsupported long double [\#1189](https://github.com/kokkos/kokkos/issues/1189)
- CUDA9: fix warning on defaulted function marking [\#1188](https://github.com/kokkos/kokkos/issues/1188)
- CUDA9: fix warnings for deprecated warp level functions [\#1187](https://github.com/kokkos/kokkos/issues/1187)
- Add CUDA 9.0 nightly testing [\#1174](https://github.com/kokkos/kokkos/issues/1174)
- {OMPI,MPICH}\_CXX hack breaks nvcc\_wrapper use case [\#1166](https://github.com/kokkos/kokkos/issues/1166)
- KOKKOS\_HAVE\_CUDA\_LAMBDA became KOKKOS\_CUDA\_USE\_LAMBDA [\#1274](https://github.com/kokkos/kokkos/issues/1274)

**Fixed bugs:**

- MinMax Reducer with tagged operator doesn't compile [\#1251](https://github.com/kokkos/kokkos/issues/1251)
- Reducers for Tagged operators give wrong answer [\#1250](https://github.com/kokkos/kokkos/issues/1250)
- Kokkos not Compatible with Big Endian Machines? [\#1235](https://github.com/kokkos/kokkos/issues/1235)
- Parallel Scan hangs forever on BG/Q [\#1234](https://github.com/kokkos/kokkos/issues/1234)
- Threads backend doesn't compile with Clang on OS X [\#1232](https://github.com/kokkos/kokkos/issues/1232)
- $\(shell date\) needs quote [\#1264](https://github.com/kokkos/kokkos/issues/1264)
- Unqualified parallel\_for call conflicts with user-defined parallel\_for [\#1219](https://github.com/kokkos/kokkos/issues/1219)
- KokkosAlgorithms: CMake issue in unit tests [\#1212](https://github.com/kokkos/kokkos/issues/1212)
- Intel 18 Error: "simd pragma has been deprecated" [\#1210](https://github.com/kokkos/kokkos/issues/1210)
- Memory leak in Kokkos::initialize [\#1194](https://github.com/kokkos/kokkos/issues/1194)
- CUDA9: compiler error with static assert template arguments [\#1190](https://github.com/kokkos/kokkos/issues/1190)
- Kokkos::Serial::is\_initialized returns always true [\#1184](https://github.com/kokkos/kokkos/issues/1184)
- Triple nested parallelism still fails on bowman [\#1093](https://github.com/kokkos/kokkos/issues/1093)
- OpenMP openmp.range on Develop Runs Forever on POWER7+ with RHEL7 and GCC4.8.5 [\#995](https://github.com/kokkos/kokkos/issues/995)
- Rendezvous performance at global scope [\#985](https://github.com/kokkos/kokkos/issues/985)


## [2.04.11](https://github.com/kokkos/kokkos/tree/2.04.11) (2017-10-28)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.04.04...2.04.11)

**Implemented enhancements:**

- Add Subview pattern. [\#648](https://github.com/kokkos/kokkos/issues/648)
- Add Kokkos "global" is\_initialized [\#1060](https://github.com/kokkos/kokkos/issues/1060)
- Add create\_mirror\_view\_and\_copy [\#1161](https://github.com/kokkos/kokkos/issues/1161)
- Add KokkosConcepts SpaceAccessibility function [\#1092](https://github.com/kokkos/kokkos/issues/1092)
- Option to Disable Initialize Warnings [\#1142](https://github.com/kokkos/kokkos/issues/1142)
- Mature task-DAG capability [\#320](https://github.com/kokkos/kokkos/issues/320)
- Promote Work DAG from experimental [\#1126](https://github.com/kokkos/kokkos/issues/1126)
- Implement new WorkGraph push/pop [\#1108](https://github.com/kokkos/kokkos/issues/1108)
- Kokkos\_ENABLE\_Cuda\_Lambda should default ON [\#1101](https://github.com/kokkos/kokkos/issues/1101)
- Add multidimensional parallel for example and improve unit test [\#1064](https://github.com/kokkos/kokkos/issues/1064)
- Fix ROCm:  Performance tests not building [\#1038](https://github.com/kokkos/kokkos/issues/1038)
- Make KOKKOS\_ALIGN\_SIZE a configure-time option [\#1004](https://github.com/kokkos/kokkos/issues/1004)
- Make alignment consistent [\#809](https://github.com/kokkos/kokkos/issues/809)
- Improve subview construction on Cuda backend [\#615](https://github.com/kokkos/kokkos/issues/615)

**Fixed bugs:**

- Kokkos::vector fixes for application [\#1134](https://github.com/kokkos/kokkos/issues/1134)
- DynamicView non-power of two value\_type [\#1177](https://github.com/kokkos/kokkos/issues/1177)
- Memory pool bug [\#1154](https://github.com/kokkos/kokkos/issues/1154)
- Cuda launch bounds performance regression bug [\#1140](https://github.com/kokkos/kokkos/issues/1140)
- Significant performance regression in LAMMPS after updating Kokkos [\#1139](https://github.com/kokkos/kokkos/issues/1139)
- CUDA compile error [\#1128](https://github.com/kokkos/kokkos/issues/1128)
- MDRangePolicy neg idx test failure in debug mode [\#1113](https://github.com/kokkos/kokkos/issues/1113)
- subview construction on Cuda backend [\#615](https://github.com/kokkos/kokkos/issues/615)

## [2.04.04](https://github.com/kokkos/kokkos/tree/2.04.04) (2017-09-11)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.04.00...2.04.04)

**Implemented enhancements:**

- OpenMP partition: set number of threads on nested level [\#1082](https://github.com/kokkos/kokkos/issues/1082)
- Add StaticCrsGraph row\(\) method [\#1071](https://github.com/kokkos/kokkos/issues/1071)
- Enhance Kokkos complex operator overloading [\#1052](https://github.com/kokkos/kokkos/issues/1052)
- Tell Trilinos packages about host+device lambda [\#1019](https://github.com/kokkos/kokkos/issues/1019)
- Function markup for defaulted class members [\#952](https://github.com/kokkos/kokkos/issues/952)
- Add deterministic random number generator [\#857](https://github.com/kokkos/kokkos/issues/857)

**Fixed bugs:**

- Fix reduction\_identity\<T\>::max for floating point numbers [\#1048](https://github.com/kokkos/kokkos/issues/1048)
- Fix MD iteration policy ignores lower bound on GPUs [\#1041](https://github.com/kokkos/kokkos/issues/1041)
- (Experimental) HBWSpace  Linking issues in KokkosKernels [\#1094](https://github.com/kokkos/kokkos/issues/1094)
- (Experimental) ROCm:  algorithms/unit\_tests test\_sort failing with segfault [\#1070](https://github.com/kokkos/kokkos/issues/1070)

## [2.04.00](https://github.com/kokkos/kokkos/tree/2.04.00) (2017-08-16)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.03.13...2.04.00)

**Implemented enhancements:**

- Added ROCm backend to support AMD GPUs
- Kokkos::complex\<T\> behaves slightly differently from std::complex\<T\> [\#1011](https://github.com/kokkos/kokkos/issues/1011)
- Kokkos::Experimental::Crs constructor arguments were in the wrong order [\#992](https://github.com/kokkos/kokkos/issues/992)
- Work graph construction ease-of-use (one lambda for count and fill) [\#991](https://github.com/kokkos/kokkos/issues/991)
- when\_all returns pointer of futures (improved interface) [\#990](https://github.com/kokkos/kokkos/issues/990)
- Allow assignment of LayoutLeft to LayoutRight or vice versa for rank-0 Views [\#594](https://github.com/kokkos/kokkos/issues/594)
- Changed the meaning of Kokkos\_ENABLE\_CXX11\_DISPATCH\_LAMBDA [\#1035](https://github.com/kokkos/kokkos/issues/1035)

**Fixed bugs:**

- memory pool default constructor does not properly set member variables. [\#1007](https://github.com/kokkos/kokkos/issues/1007)

## [2.03.13](https://github.com/kokkos/kokkos/tree/2.03.13) (2017-07-27)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.03.05...2.03.13)

**Implemented enhancements:**

- Disallow enabling both OpenMP and Threads in the same executable [\#406](https://github.com/kokkos/kokkos/issues/406)
- Make Kokkos::OpenMP respect OMP environment even if hwloc is available [\#630](https://github.com/kokkos/kokkos/issues/630)
- Improve Atomics Performance on KNL/Broadwell where PREFETCHW/RFO is Available [\#898](https://github.com/kokkos/kokkos/issues/898)
- Kokkos::resize should test whether dimensions have changed before resizing [\#904](https://github.com/kokkos/kokkos/issues/904)
- Develop performance-regression/acceptance tests [\#737](https://github.com/kokkos/kokkos/issues/737)
- Make the deep\_copy Profiling hook a start/end system [\#890](https://github.com/kokkos/kokkos/issues/890)
- Add deep\_copy Profiling hook [\#843](https://github.com/kokkos/kokkos/issues/843)
- Append tag name to parallel construct name for Profiling [\#842](https://github.com/kokkos/kokkos/issues/842)
- Add view label to `View bounds error` message for CUDA backend [\#870](https://github.com/kokkos/kokkos/issues/870)
- Disable printing the loaded profiling library [\#824](https://github.com/kokkos/kokkos/issues/824)
- "Declared but never referenced" warnings [\#853](https://github.com/kokkos/kokkos/issues/853)
- Warnings about lock\_address\_cuda\_space [\#852](https://github.com/kokkos/kokkos/issues/852)
- WorkGraph execution policy [\#771](https://github.com/kokkos/kokkos/issues/771)
- Simplify makefiles by guarding compilation with appropriate KOKKOS\_ENABLE\_\#\#\# macros [\#716](https://github.com/kokkos/kokkos/issues/716)
- Cmake build: wrong include install directory [\#668](https://github.com/kokkos/kokkos/issues/668)
- Derived View type and allocation [\#566](https://github.com/kokkos/kokkos/issues/566)
- Fix Compiler warnings when compiling core unit tests for Cuda [\#214](https://github.com/kokkos/kokkos/issues/214)

**Fixed bugs:**

- Out-of-bounds read in Kokkos\_Layout.hpp [\#975](https://github.com/kokkos/kokkos/issues/975)
- CudaClang: Fix failing test with Clang 4.0 [\#941](https://github.com/kokkos/kokkos/issues/941)
- Respawn when memory pool allocation fails \(not available memory\) [\#940](https://github.com/kokkos/kokkos/issues/940)
- Memory pool aborts on zero allocation request, returns NULL for \< minimum [\#939](https://github.com/kokkos/kokkos/issues/939)
- Error with TaskScheduler query of underlying memory pool [\#917](https://github.com/kokkos/kokkos/issues/917)
- Profiling::\*Callee static variables declared in header [\#863](https://github.com/kokkos/kokkos/issues/863)
- calling \*Space::name\(\) causes compile error [\#862](https://github.com/kokkos/kokkos/issues/862)
- bug in Profiling::deallocateData [\#860](https://github.com/kokkos/kokkos/issues/860)
- task\_depend test failing, CUDA 8.0 + Pascal + RDC [\#829](https://github.com/kokkos/kokkos/issues/829)
- \[develop branch\] Standalone cmake issues [\#826](https://github.com/kokkos/kokkos/issues/826)
- Kokkos CUDA failes to compile with OMPI\_CXX and MPICH\_CXX wrappers [\#776](https://github.com/kokkos/kokkos/issues/776)
- Task Team reduction on Pascal [\#767](https://github.com/kokkos/kokkos/issues/767)
- CUDA stack overflow with TaskDAG test [\#758](https://github.com/kokkos/kokkos/issues/758)
- TeamVector test on Cuda [\#670](https://github.com/kokkos/kokkos/issues/670)
- Clang 4.0 Cuda Build broken again [\#560](https://github.com/kokkos/kokkos/issues/560)


## [2.03.05](https://github.com/kokkos/kokkos/tree/2.03.05) (2017-05-27)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.03.00...2.03.05)

**Implemented enhancements:**

- Harmonize Custom Reductions over nesting levels [\#802](https://github.com/kokkos/kokkos/issues/802)
- Prevent users directly including KokkosCore\_config.h [\#815](https://github.com/kokkos/kokkos/issues/815)
- DualView aborts on concurrent host/device modify \(in debug mode\) [\#814](https://github.com/kokkos/kokkos/issues/814)
- Abort when running on a NVIDIA CC5.0 or higher architecture with code compiled for CC \< 5.0 [\#813](https://github.com/kokkos/kokkos/issues/813)
- Add "name" function to ExecSpaces [\#806](https://github.com/kokkos/kokkos/issues/806)
- Allow null Future in task spawn dependences [\#795](https://github.com/kokkos/kokkos/issues/795)
- Add Unit Tests for Kokkos::complex [\#785](https://github.com/kokkos/kokkos/issues/785)
- Add pow function for Kokkos::complex [\#784](https://github.com/kokkos/kokkos/issues/784)
- Square root of a complex [\#729](https://github.com/kokkos/kokkos/issues/729)
- Command line processing of --threads argument prevents users from having any commandline arguments starting with --threads [\#760](https://github.com/kokkos/kokkos/issues/760)
- Protected deprecated API with appropriate macro [\#756](https://github.com/kokkos/kokkos/issues/756)
- Allow task scheduler memory pool to be used by tasks [\#747](https://github.com/kokkos/kokkos/issues/747)
- View bounds checking on host-side performance: constructing a std::string [\#723](https://github.com/kokkos/kokkos/issues/723)
- Add check for AppleClang as compiler distinct from check for Clang. [\#705](https://github.com/kokkos/kokkos/issues/705)
- Uninclude source files for specific configurations to prevent link warning. [\#701](https://github.com/kokkos/kokkos/issues/701)
- Add --small option to snapshot script [\#697](https://github.com/kokkos/kokkos/issues/697)
- CMake Standalone Support [\#674](https://github.com/kokkos/kokkos/issues/674)
- CMake build unit test and install [\#808](https://github.com/kokkos/kokkos/issues/808)
- CMake: Fix having kokkos as a subdirectory in a pure cmake project [\#629](https://github.com/kokkos/kokkos/issues/629)
- Tribits macro assumes build directory is in top level source directory [\#654](https://github.com/kokkos/kokkos/issues/654)
- Use bin/nvcc\_wrapper, not config/nvcc\_wrapper [\#562](https://github.com/kokkos/kokkos/issues/562)
- Allow MemoryPool::allocate\(\) to be called from multiple threads per warp. [\#487](https://github.com/kokkos/kokkos/issues/487)
- Allow MemoryPool::allocate\\(\\) to be called from multiple threads per warp. [\#487](https://github.com/kokkos/kokkos/issues/487)
- Move OpenMP 4.5 OpenMPTarget backend into Develop [\#456](https://github.com/kokkos/kokkos/issues/456)
- Testing on ARM testbed [\#288](https://github.com/kokkos/kokkos/issues/288)

**Fixed bugs:**

- Fix label in OpenMP parallel\_reduce verify\_initialized [\#834](https://github.com/kokkos/kokkos/issues/834)
- TeamScratch Level 1 on Cuda hangs [\#820](https://github.com/kokkos/kokkos/issues/820)
- \[bug\] memory pool. [\#786](https://github.com/kokkos/kokkos/issues/786)
- Some Reduction Tests fail on Intel 18 with aggressive vectorization on [\#774](https://github.com/kokkos/kokkos/issues/774)
- Error copying dynamic view on copy of memory pool [\#773](https://github.com/kokkos/kokkos/issues/773)
- CUDA stack overflow with TaskDAG test [\#758](https://github.com/kokkos/kokkos/issues/758)
- ThreadVectorRange Customized Reduction Bug [\#739](https://github.com/kokkos/kokkos/issues/739)
- set\_scratch\_size overflows  [\#726](https://github.com/kokkos/kokkos/issues/726)
- Get wrong results for compiler checks in Makefile on OS X. [\#706](https://github.com/kokkos/kokkos/issues/706)
- Fix check if multiple host architectures enabled. [\#702](https://github.com/kokkos/kokkos/issues/702)
- Threads Backend Does not Pass on Cray Compilers [\#609](https://github.com/kokkos/kokkos/issues/609)
- Rare bug in memory pool where allocation can finish on superblock in empty state [\#452](https://github.com/kokkos/kokkos/issues/452)
- LDFLAGS in core/unit\_test/Makefile: potential "undefined reference" to pthread lib [\#148](https://github.com/kokkos/kokkos/issues/148)

## [2.03.00](https://github.com/kokkos/kokkos/tree/2.03.00) (2017-04-25)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.02.15...2.03.00)

**Implemented enhancements:**

- UnorderedMap: make it accept Devices or MemorySpaces [\#711](https://github.com/kokkos/kokkos/issues/711)
- sort to accept DynamicView and \[begin,end\) indices [\#691](https://github.com/kokkos/kokkos/issues/691)
- ENABLE Macros should only be used via \#ifdef or \#if defined [\#675](https://github.com/kokkos/kokkos/issues/675)
- Remove impl/Kokkos\_Synchronic\_\* [\#666](https://github.com/kokkos/kokkos/issues/666)
- Turning off IVDEP for Intel 14.  [\#638](https://github.com/kokkos/kokkos/issues/638)
- Using an installed Kokkos in a target application using CMake [\#633](https://github.com/kokkos/kokkos/issues/633)
- Create Kokkos Bill of Materials [\#632](https://github.com/kokkos/kokkos/issues/632)
- MDRangePolicy and tagged evaluators [\#547](https://github.com/kokkos/kokkos/issues/547)
- Add PGI support [\#289](https://github.com/kokkos/kokkos/issues/289)

**Fixed bugs:**

- Output from PerTeam fails [\#733](https://github.com/kokkos/kokkos/issues/733)
- Cuda: architecture flag not added to link line [\#688](https://github.com/kokkos/kokkos/issues/688)
- Getting large chunks of memory for a thread team in a universal way [\#664](https://github.com/kokkos/kokkos/issues/664)
- Kokkos RNG normal\(\) function hangs for small seed value [\#655](https://github.com/kokkos/kokkos/issues/655)
- Kokkos Tests Errors on Shepard/HSW Builds [\#644](https://github.com/kokkos/kokkos/issues/644)

## [2.02.15](https://github.com/kokkos/kokkos/tree/2.02.15) (2017-02-10)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.02.07...2.02.15)

**Implemented enhancements:**

- Containers: Adding block partitioning to StaticCrsGraph [\#625](https://github.com/kokkos/kokkos/issues/625)
- Kokkos Make System can induce Errors on Cray Volta System [\#610](https://github.com/kokkos/kokkos/issues/610)
- OpenMP: error out if KOKKOS\_HAVE\_OPENMP is defined but not \_OPENMP [\#605](https://github.com/kokkos/kokkos/issues/605)
- CMake: fix standalone build with tests [\#604](https://github.com/kokkos/kokkos/issues/604)
- Change README \(that GitHub shows when opening Kokkos project page\) to tell users how to submit PRs [\#597](https://github.com/kokkos/kokkos/issues/597)
- Add correctness testing for all operators of Atomic View [\#420](https://github.com/kokkos/kokkos/issues/420)
- Allow assignment of Views with compatible memory spaces [\#290](https://github.com/kokkos/kokkos/issues/290)
- Build only one version of Kokkos library for tests [\#213](https://github.com/kokkos/kokkos/issues/213)
- Clean out old KOKKOS\_HAVE\_CXX11 macros clauses [\#156](https://github.com/kokkos/kokkos/issues/156)
- Harmonize Macro names [\#150](https://github.com/kokkos/kokkos/issues/150)

**Fixed bugs:**

- Cray and PGI: Kokkos\_Parallel\_Reduce [\#634](https://github.com/kokkos/kokkos/issues/634)
- Kokkos Make System can induce Errors on Cray Volta System [\#610](https://github.com/kokkos/kokkos/issues/610)
- Normal\(\) function random number generator doesn't give the expected distribution [\#592](https://github.com/kokkos/kokkos/issues/592)

## [2.02.07](https://github.com/kokkos/kokkos/tree/2.02.07) (2016-12-16)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.02.01...2.02.07)

**Implemented enhancements:**

- Add CMake option to enable Cuda Lambda support [\#589](https://github.com/kokkos/kokkos/issues/589)
- Add CMake option to enable Cuda RDC support [\#588](https://github.com/kokkos/kokkos/issues/588)
- Add Initial Intel Sky Lake Xeon-HPC Compiler Support to Kokkos Make System [\#584](https://github.com/kokkos/kokkos/issues/584)
- Building Tutorial Examples  [\#582](https://github.com/kokkos/kokkos/issues/582)
- Internal way for using ThreadVectorRange without TeamHandle  [\#574](https://github.com/kokkos/kokkos/issues/574)
- Testing: Add testing for uvm and rdc [\#571](https://github.com/kokkos/kokkos/issues/571)
- Profiling: Add Memory Tracing and Region Markers [\#557](https://github.com/kokkos/kokkos/issues/557)
- nvcc\_wrapper not installed with Kokkos built with CUDA through CMake [\#543](https://github.com/kokkos/kokkos/issues/543)
- Improve DynRankView debug check [\#541](https://github.com/kokkos/kokkos/issues/541)
- Benchmarks: Add Gather benchmark [\#536](https://github.com/kokkos/kokkos/issues/536)
- Testing: add spot\_check option to test\_all\_sandia [\#535](https://github.com/kokkos/kokkos/issues/535)
- Deprecate Kokkos::Impl::VerifyExecutionCanAccessMemorySpace [\#527](https://github.com/kokkos/kokkos/issues/527)
- Add AtomicAdd support for 64bit float for Pascal [\#522](https://github.com/kokkos/kokkos/issues/522)
- Add Restrict and Aligned memory trait [\#517](https://github.com/kokkos/kokkos/issues/517)
- Kokkos Tests are Not Run using Compiler Optimization [\#501](https://github.com/kokkos/kokkos/issues/501)
- Add support for clang 3.7 w/ openmp backend [\#393](https://github.com/kokkos/kokkos/issues/393)
- Provide an error throw class [\#79](https://github.com/kokkos/kokkos/issues/79)

**Fixed bugs:**

- Cuda UVM Allocation test broken with UVM as default space [\#586](https://github.com/kokkos/kokkos/issues/586)
- Bug \(develop branch only\): multiple tests are now failing when forcing uvm usage. [\#570](https://github.com/kokkos/kokkos/issues/570)
- Error in generate\_makefile.sh for Kokkos when Compiler is Empty String/Fails [\#568](https://github.com/kokkos/kokkos/issues/568)
- XL 13.1.4 incorrect C++11 flag [\#553](https://github.com/kokkos/kokkos/issues/553)
- Improve DynRankView debug check [\#541](https://github.com/kokkos/kokkos/issues/541)
- Installing Library on MAC broken due to cp -u [\#539](https://github.com/kokkos/kokkos/issues/539)
- Intel Nightly Testing with Debug enabled fails [\#534](https://github.com/kokkos/kokkos/issues/534)

## [2.02.01](https://github.com/kokkos/kokkos/tree/2.02.01) (2016-11-01)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.02.00...2.02.01)

**Implemented enhancements:**

- Add Changelog generation to our process. [\#506](https://github.com/kokkos/kokkos/issues/506)

**Fixed bugs:**

- Test scratch\_request fails in Serial with Debug enabled [\#520](https://github.com/kokkos/kokkos/issues/520)
- Bug In BoundsCheck for DynRankView [\#516](https://github.com/kokkos/kokkos/issues/516)

## [2.02.00](https://github.com/kokkos/kokkos/tree/2.02.00) (2016-10-30)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.01.10...2.02.00)

**Implemented enhancements:**

- Add PowerPC assembly for grabbing clock register in memory pool [\#511](https://github.com/kokkos/kokkos/issues/511)
- Add GCC 6.x support [\#508](https://github.com/kokkos/kokkos/issues/508)
- Test install and build against installed library [\#498](https://github.com/kokkos/kokkos/issues/498)
- Makefile.kokkos adds expt-extended-lambda to cuda build with clang [\#490](https://github.com/kokkos/kokkos/issues/490)
- Add top-level makefile option to just test kokkos-core unit-test [\#485](https://github.com/kokkos/kokkos/issues/485)
- Split and harmonize Object Files of Core UnitTests to increase build parallelism [\#484](https://github.com/kokkos/kokkos/issues/484)
- LayoutLeft to LayoutLeft subview for 3D and 4D views [\#473](https://github.com/kokkos/kokkos/issues/473)
- Add official Cuda 8.0 support [\#468](https://github.com/kokkos/kokkos/issues/468)
- Allow C++1Z Flag for Class Lambda capture [\#465](https://github.com/kokkos/kokkos/issues/465)
- Add Clang 4.0+ compilation of Cuda code [\#455](https://github.com/kokkos/kokkos/issues/455)
- Possible Issue with Intel 17.0.098 and GCC 6.1.0 in Develop Branch [\#445](https://github.com/kokkos/kokkos/issues/445)
- Add name of view to "View bounds error" [\#432](https://github.com/kokkos/kokkos/issues/432)
- Move Sort Binning Operators into Kokkos namespace [\#421](https://github.com/kokkos/kokkos/issues/421)
- TaskPolicy - generate error when attempt to use uninitialized  [\#396](https://github.com/kokkos/kokkos/issues/396)
- Import WithoutInitializing and AllowPadding into Kokkos namespace [\#325](https://github.com/kokkos/kokkos/issues/325)
- TeamThreadRange requires begin, end to be the same type [\#305](https://github.com/kokkos/kokkos/issues/305)
- CudaUVMSpace should track \# allocations, due to CUDA limit on \# UVM allocations [\#300](https://github.com/kokkos/kokkos/issues/300)
- Remove old View and its infrastructure [\#259](https://github.com/kokkos/kokkos/issues/259)

**Fixed bugs:**

- Bug in TestCuda\_Other.cpp: most likely assembly inserted into Device code [\#515](https://github.com/kokkos/kokkos/issues/515)
- Cuda Compute Capability check of GPU is outdated [\#509](https://github.com/kokkos/kokkos/issues/509)
- multi\_scratch test with hwloc and pthreads seg-faults.  [\#504](https://github.com/kokkos/kokkos/issues/504)
- generate\_makefile.bash: "make install" is broken [\#503](https://github.com/kokkos/kokkos/issues/503)
- make clean in Out of Source Build/Tests Does Not Work Correctly [\#502](https://github.com/kokkos/kokkos/issues/502)
- Makefiles for test and examples have issues in Cuda when CXX is not explicitly specified [\#497](https://github.com/kokkos/kokkos/issues/497)
- Dispatch lambda test directly inside GTEST macro doesn't work with nvcc [\#491](https://github.com/kokkos/kokkos/issues/491)
- UnitTests with HWLOC enabled fail if run with mpirun bound to a single core [\#489](https://github.com/kokkos/kokkos/issues/489)
- Failing Reducer Test on Mac with Pthreads [\#479](https://github.com/kokkos/kokkos/issues/479)
- make test Dumps Error with Clang Not Found [\#471](https://github.com/kokkos/kokkos/issues/471)
- OpenMP TeamPolicy member broadcast not using correct volatile shared variable [\#424](https://github.com/kokkos/kokkos/issues/424)
- TaskPolicy - generate error when attempt to use uninitialized  [\#396](https://github.com/kokkos/kokkos/issues/396)
- New task policy implementation is pulling in old experimental code. [\#372](https://github.com/kokkos/kokkos/issues/372)
- MemoryPool unit test hangs on Power8 with GCC 6.1.0 [\#298](https://github.com/kokkos/kokkos/issues/298)

## [2.01.10](https://github.com/kokkos/kokkos/tree/2.01.10) (2016-09-27)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.01.06...2.01.10)

**Implemented enhancements:**

- Enable Profiling by default in Tribits build [\#438](https://github.com/kokkos/kokkos/issues/438)
- parallel\_reduce\(0\), parallel\_scan\(0\) unit tests [\#436](https://github.com/kokkos/kokkos/issues/436)
- data\(\)==NULL after realloc with LayoutStride [\#351](https://github.com/kokkos/kokkos/issues/351)
- Fix tutorials to track new Kokkos::View [\#323](https://github.com/kokkos/kokkos/issues/323)
- Rename team policy set\_scratch\_size. [\#195](https://github.com/kokkos/kokkos/issues/195)

**Fixed bugs:**

- Possible Issue with Intel 17.0.098 and GCC 6.1.0 in Develop Branch [\#445](https://github.com/kokkos/kokkos/issues/445)
- Makefile spits syntax error [\#435](https://github.com/kokkos/kokkos/issues/435)
- Kokkos::sort fails for view with all the same values [\#422](https://github.com/kokkos/kokkos/issues/422)
- Generic Reducers: can't accept inline constructed reducer [\#404](https://github.com/kokkos/kokkos/issues/404)
- data\\(\\)==NULL after realloc with LayoutStride [\#351](https://github.com/kokkos/kokkos/issues/351)
- const subview of const view with compile time dimensions on Cuda backend [\#310](https://github.com/kokkos/kokkos/issues/310)
- Kokkos \(in Trilinos\) Causes Internal Compiler Error on CUDA 8.0.21-EA on POWER8 [\#307](https://github.com/kokkos/kokkos/issues/307)
- Core Oversubscription Detection Broken? [\#159](https://github.com/kokkos/kokkos/issues/159)


## [2.01.06](https://github.com/kokkos/kokkos/tree/2.01.06) (2016-09-02)
[Full Changelog](https://github.com/kokkos/kokkos/compare/2.01.00...2.01.06)

**Implemented enhancements:**

- Add "standard" reducers for lambda-supportable customized reduce [\#411](https://github.com/kokkos/kokkos/issues/411)
- TaskPolicy - single thread back-end execution [\#390](https://github.com/kokkos/kokkos/issues/390)
- Kokkos master clone tag [\#387](https://github.com/kokkos/kokkos/issues/387)
- Query memory requirements from task policy [\#378](https://github.com/kokkos/kokkos/issues/378)
- Output order of test\_atomic.cpp is confusing [\#373](https://github.com/kokkos/kokkos/issues/373)
- Missing testing for atomics [\#341](https://github.com/kokkos/kokkos/issues/341)
- Feature request for Kokkos to provide Kokkos::atomic\_fetch\_max and atomic\_fetch\_min [\#336](https://github.com/kokkos/kokkos/issues/336)
- TaskPolicy\<Cuda\> performance requires teams mapped to warps [\#218](https://github.com/kokkos/kokkos/issues/218)

**Fixed bugs:**

- Reduce with Teams broken for custom initialize [\#407](https://github.com/kokkos/kokkos/issues/407)
- Failing Kokkos build on Debian [\#402](https://github.com/kokkos/kokkos/issues/402)
- Failing Tests on NVIDIA Pascal GPUs [\#398](https://github.com/kokkos/kokkos/issues/398)
- Algorithms: fill\_random assumes dimensions fit in unsigned int [\#389](https://github.com/kokkos/kokkos/issues/389)
- Kokkos::subview with RandomAccess Memory Trait [\#385](https://github.com/kokkos/kokkos/issues/385)
- Build warning \(signed / unsigned comparison\) in Cuda implementation [\#365](https://github.com/kokkos/kokkos/issues/365)
- wrong results for a parallel\_reduce with CUDA8 / Maxwell50 [\#352](https://github.com/kokkos/kokkos/issues/352)
- Hierarchical parallelism - 3 level unit test [\#344](https://github.com/kokkos/kokkos/issues/344)
- Can I allocate a View w/ both WithoutInitializing & AllowPadding? [\#324](https://github.com/kokkos/kokkos/issues/324)
- subview View layout determination [\#309](https://github.com/kokkos/kokkos/issues/309)
- Unit tests with Cuda - Maxwell [\#196](https://github.com/kokkos/kokkos/issues/196)

## [2.01.00](https://github.com/kokkos/kokkos/tree/2.01.00) (2016-07-21)
[Full Changelog](https://github.com/kokkos/kokkos/compare/End_C++98...2.01.00)

**Implemented enhancements:**

- Edit ViewMapping so assigning Views with the same custom layout compiles when const casting [\#327](https://github.com/kokkos/kokkos/issues/327)
- DynRankView: Performance improvement for operator\(\) [\#321](https://github.com/kokkos/kokkos/issues/321)
- Interoperability between static and dynamic rank views [\#295](https://github.com/kokkos/kokkos/issues/295)
- subview member function ? [\#280](https://github.com/kokkos/kokkos/issues/280)
- Inter-operatibility between View and DynRankView. [\#245](https://github.com/kokkos/kokkos/issues/245)
- \(Trilinos\) build warning in atomic\_assign, with Kokkos::complex [\#177](https://github.com/kokkos/kokkos/issues/177)
- View\<\>::shmem\_size should runtime check for number of arguments equal to rank [\#176](https://github.com/kokkos/kokkos/issues/176)
- Custom reduction join via lambda argument [\#99](https://github.com/kokkos/kokkos/issues/99)
- DynRankView with 0 dimensions passed in at construction [\#293](https://github.com/kokkos/kokkos/issues/293)
- Inject view\_alloc and friends into Kokkos namespace [\#292](https://github.com/kokkos/kokkos/issues/292)
- Less restrictive TeamPolicy reduction on Cuda [\#286](https://github.com/kokkos/kokkos/issues/286)
- deep\_copy using remap with source execution space [\#267](https://github.com/kokkos/kokkos/issues/267)
- Suggestion:  Enable opt-in L1 caching via nvcc-wrapper [\#261](https://github.com/kokkos/kokkos/issues/261)
- More flexible create\_mirror functions [\#260](https://github.com/kokkos/kokkos/issues/260)
- Rename View::memory\_span to View::required\_allocation\_size [\#256](https://github.com/kokkos/kokkos/issues/256)
- Use of subviews and views with compile-time dimensions [\#237](https://github.com/kokkos/kokkos/issues/237)
- Use of subviews and views with compile-time dimensions [\#237](https://github.com/kokkos/kokkos/issues/237)
- Kokkos::Timer [\#234](https://github.com/kokkos/kokkos/issues/234)
- Fence CudaUVMSpace allocations [\#230](https://github.com/kokkos/kokkos/issues/230)
- View::operator\(\) accept std::is\_integral and std::is\_enum [\#227](https://github.com/kokkos/kokkos/issues/227)
- Allocating zero size View [\#216](https://github.com/kokkos/kokkos/issues/216)
- Thread scalable memory pool [\#212](https://github.com/kokkos/kokkos/issues/212)
- Add a way to disable memory leak output [\#194](https://github.com/kokkos/kokkos/issues/194)
- Kokkos exec space init should init Kokkos profiling [\#192](https://github.com/kokkos/kokkos/issues/192)
- Runtime rank wrapper for View [\#189](https://github.com/kokkos/kokkos/issues/189)
- Profiling Interface [\#158](https://github.com/kokkos/kokkos/issues/158)
- Fix View assignment \(of managed to unmanaged\) [\#153](https://github.com/kokkos/kokkos/issues/153)
- Add unit test for assignment of managed View to unmanaged View [\#152](https://github.com/kokkos/kokkos/issues/152)
- Check for oversubscription of threads with MPI in Kokkos::initialize [\#149](https://github.com/kokkos/kokkos/issues/149)
- Dynamic resizeable 1dimensional view [\#143](https://github.com/kokkos/kokkos/issues/143)
- Develop TaskPolicy for CUDA [\#142](https://github.com/kokkos/kokkos/issues/142)
- New View : Test Compilation Downstream [\#138](https://github.com/kokkos/kokkos/issues/138)
- New View Implementation [\#135](https://github.com/kokkos/kokkos/issues/135)
- Add variant of subview that lets users add traits [\#134](https://github.com/kokkos/kokkos/issues/134)
- NVCC-WRAPPER: Add --host-only flag [\#121](https://github.com/kokkos/kokkos/issues/121)
- Address gtest issue with TriBITS Kokkos build outside of Trilinos [\#117](https://github.com/kokkos/kokkos/issues/117)
- Make tests pass with -expt-extended-lambda on CUDA [\#108](https://github.com/kokkos/kokkos/issues/108)
- Dynamic scheduling for parallel\_for and parallel\_reduce [\#106](https://github.com/kokkos/kokkos/issues/106)
- Runtime or compile time error when reduce functor's join is not properly specified as const member function or with volatile arguments [\#105](https://github.com/kokkos/kokkos/issues/105)
- Error out when the number of threads is modified after kokkos is initialized [\#104](https://github.com/kokkos/kokkos/issues/104)
- Porting to POWER and remove assumption of X86 default [\#103](https://github.com/kokkos/kokkos/issues/103)
- Dynamic scheduling option for RangePolicy [\#100](https://github.com/kokkos/kokkos/issues/100)
- SharedMemory Support for Lambdas [\#81](https://github.com/kokkos/kokkos/issues/81)
- Recommended TeamSize for Lambdas [\#80](https://github.com/kokkos/kokkos/issues/80)
- Add Aggressive Vectorization Compilation mode [\#72](https://github.com/kokkos/kokkos/issues/72)
- Dynamic scheduling team execution policy [\#53](https://github.com/kokkos/kokkos/issues/53)
- UVM allocations in multi-GPU systems [\#50](https://github.com/kokkos/kokkos/issues/50)
- Synchronic in Kokkos::Impl [\#44](https://github.com/kokkos/kokkos/issues/44)
- index and dimension types in for loops [\#28](https://github.com/kokkos/kokkos/issues/28)
- Subview assign of 1D Strided with stride 1 to LayoutLeft/Right [\#1](https://github.com/kokkos/kokkos/issues/1)

**Fixed bugs:**

- misspelled variable name in Kokkos\_Atomic\_Fetch + missing unit tests [\#340](https://github.com/kokkos/kokkos/issues/340)
- seg fault Kokkos::Impl::CudaInternal::print\_configuration [\#338](https://github.com/kokkos/kokkos/issues/338)
- Clang compiler error with named parallel\_reduce, tags, and TeamPolicy. [\#335](https://github.com/kokkos/kokkos/issues/335)
- Shared Memory Allocation Error at parallel\_reduce [\#311](https://github.com/kokkos/kokkos/issues/311)
- DynRankView: Fix resize and realloc [\#303](https://github.com/kokkos/kokkos/issues/303)
- Scratch memory and dynamic scheduling [\#279](https://github.com/kokkos/kokkos/issues/279)
- MemoryPool infinite loop when out of memory [\#312](https://github.com/kokkos/kokkos/issues/312)
- Kokkos DynRankView changes break Sacado and Panzer [\#299](https://github.com/kokkos/kokkos/issues/299)
- MemoryPool fails to compile on non-cuda non-x86 [\#297](https://github.com/kokkos/kokkos/issues/297)
- Random Number Generator Fix [\#296](https://github.com/kokkos/kokkos/issues/296)
- View template parameter ordering Bug [\#282](https://github.com/kokkos/kokkos/issues/282)
- Serial task policy broken. [\#281](https://github.com/kokkos/kokkos/issues/281)
- deep\_copy with LayoutStride should not memcpy [\#262](https://github.com/kokkos/kokkos/issues/262)
- DualView::need\_sync should be a const method [\#248](https://github.com/kokkos/kokkos/issues/248)
- Arbitrary-sized atomics on GPUs broken; loop forever [\#238](https://github.com/kokkos/kokkos/issues/238)
- boolean reduction value\_type changes answer [\#225](https://github.com/kokkos/kokkos/issues/225)
- Custom init\(\) function for parallel\_reduce with array value\_type [\#210](https://github.com/kokkos/kokkos/issues/210)
- unit\_test Makefile is Broken - Recursively Calls itself until Machine Apocalypse. [\#202](https://github.com/kokkos/kokkos/issues/202)
- nvcc\_wrapper Does Not Support  -Xcompiler \<compiler option\> [\#198](https://github.com/kokkos/kokkos/issues/198)
- Kokkos exec space init should init Kokkos profiling [\#192](https://github.com/kokkos/kokkos/issues/192)
- Kokkos Threads Backend impl\_shared\_alloc Broken on Intel 16.1 \(Shepard Haswell\) [\#186](https://github.com/kokkos/kokkos/issues/186)
- pthread back end hangs if used uninitialized [\#182](https://github.com/kokkos/kokkos/issues/182)
- parallel\_reduce of size 0, not calling init/join [\#175](https://github.com/kokkos/kokkos/issues/175)
- Bug in Threads with OpenMP enabled [\#173](https://github.com/kokkos/kokkos/issues/173)
- KokkosExp\_SharedAlloc, m\_team\_work\_index inaccessible [\#166](https://github.com/kokkos/kokkos/issues/166)
- 128-bit CAS without Assembly Broken? [\#161](https://github.com/kokkos/kokkos/issues/161)
- fatal error: Cuda/Kokkos\_Cuda\_abort.hpp: No such file or directory [\#157](https://github.com/kokkos/kokkos/issues/157)
- Power8: Fix OpenMP backend [\#139](https://github.com/kokkos/kokkos/issues/139)
- Data race in Kokkos OpenMP initialization [\#131](https://github.com/kokkos/kokkos/issues/131)
- parallel\_launch\_local\_memory and cuda 7.5 [\#125](https://github.com/kokkos/kokkos/issues/125)
- Resize can fail with Cuda due to asynchronous dispatch [\#119](https://github.com/kokkos/kokkos/issues/119)
- Qthread taskpolicy initialization bug. [\#92](https://github.com/kokkos/kokkos/issues/92)
- Windows: sys/mman.h [\#89](https://github.com/kokkos/kokkos/issues/89)
- Windows: atomic\_fetch\_sub\(\) [\#88](https://github.com/kokkos/kokkos/issues/88)
- Windows: snprintf [\#87](https://github.com/kokkos/kokkos/issues/87)
- Parallel\_Reduce with TeamPolicy and league size of 0 returns garbage [\#85](https://github.com/kokkos/kokkos/issues/85)
- Throw with Cuda when using \(2D\) team\_policy parallel\_reduce with less than a warp size [\#76](https://github.com/kokkos/kokkos/issues/76)
- Scalar views don't work with Kokkos::Atomic memory trait [\#69](https://github.com/kokkos/kokkos/issues/69)
- Reduce the number of threads per team for Cuda [\#63](https://github.com/kokkos/kokkos/issues/63)
- Named Kernels fail for reductions with CUDA [\#60](https://github.com/kokkos/kokkos/issues/60)
- Kokkos View dimension\_\(\) for long returning unsigned int [\#20](https://github.com/kokkos/kokkos/issues/20)
- atomic test hangs with LLVM [\#6](https://github.com/kokkos/kokkos/issues/6)
- OpenMP Test should set omp\_set\_num\_threads to 1 [\#4](https://github.com/kokkos/kokkos/issues/4)

**Closed issues:**

- develop branch broken with CUDA 8 and --expt-extended-lambda  [\#354](https://github.com/kokkos/kokkos/issues/354)
- --arch=KNL with Intel 2016 build failure [\#349](https://github.com/kokkos/kokkos/issues/349)
- Error building with Cuda when passing -DKOKKOS\_CUDA\_USE\_LAMBDA to generate\_makefile.bash [\#343](https://github.com/kokkos/kokkos/issues/343)
- Can I safely use int indices in a 2-D View with capacity \> 2B? [\#318](https://github.com/kokkos/kokkos/issues/318)
- Kokkos::ViewAllocateWithoutInitializing is not working [\#317](https://github.com/kokkos/kokkos/issues/317)
- Intel build on Mac OS X [\#277](https://github.com/kokkos/kokkos/issues/277)
- deleted [\#271](https://github.com/kokkos/kokkos/issues/271)
- Broken Mira build [\#268](https://github.com/kokkos/kokkos/issues/268)
- 32-bit build [\#246](https://github.com/kokkos/kokkos/issues/246)
- parallel\_reduce with RDC crashes linker [\#232](https://github.com/kokkos/kokkos/issues/232)
- build of Kokkos\_Sparse\_MV\_impl\_spmv\_Serial.cpp.o fails if you use nvcc and have cuda disabled [\#209](https://github.com/kokkos/kokkos/issues/209)
- Kokkos Serial execution space is not tested with TeamPolicy. [\#207](https://github.com/kokkos/kokkos/issues/207)
- Unit test failure on Hansen  KokkosCore\_UnitTest\_Cuda\_MPI\_1 [\#200](https://github.com/kokkos/kokkos/issues/200)
- nvcc compiler warning: calling a \_\_host\_\_ function from a \_\_host\_\_ \_\_device\_\_ function is not allowed [\#180](https://github.com/kokkos/kokkos/issues/180)
- Intel 15 build error with defaulted "move" operators [\#171](https://github.com/kokkos/kokkos/issues/171)
- missing libkokkos.a during Trilinos 12.4.2 build, yet other libkokkos\*.a libs are there [\#165](https://github.com/kokkos/kokkos/issues/165)
- Tie atomic updates to execution space or even to thread team? \(speculation\) [\#144](https://github.com/kokkos/kokkos/issues/144)
- New View: Compiletime/size Test [\#137](https://github.com/kokkos/kokkos/issues/137)
- New View : Performance Test [\#136](https://github.com/kokkos/kokkos/issues/136)
- Signed/unsigned  comparison warning in CUDA parallel [\#130](https://github.com/kokkos/kokkos/issues/130)
- Kokkos::complex: Need op\* w/ std::complex & real [\#126](https://github.com/kokkos/kokkos/issues/126)
- Use uintptr\_t for casting pointers [\#110](https://github.com/kokkos/kokkos/issues/110)
- Default thread mapping behavior between P and Q threads. [\#91](https://github.com/kokkos/kokkos/issues/91)
- Windows: Atomic\_Fetch\_Exchange\(\) return type [\#90](https://github.com/kokkos/kokkos/issues/90)
- Synchronic unit test is way too long [\#84](https://github.com/kokkos/kokkos/issues/84)
- nvcc\_wrapper -\> $\(NVCC\_WRAPPER\) [\#42](https://github.com/kokkos/kokkos/issues/42)
- Check compiler version and print helpful message [\#39](https://github.com/kokkos/kokkos/issues/39)
- Kokkos shared memory on Cuda uses a lot of registers [\#31](https://github.com/kokkos/kokkos/issues/31)
- Can not pass unit test `cuda.space` without a GT 720 [\#25](https://github.com/kokkos/kokkos/issues/25)
- Makefile.kokkos lacks bounds checking option that CMake has [\#24](https://github.com/kokkos/kokkos/issues/24)
- Kokkos can not complete unit tests with CUDA UVM enabled [\#23](https://github.com/kokkos/kokkos/issues/23)
- Simplify teams + shared memory histogram example to remove vectorization [\#21](https://github.com/kokkos/kokkos/issues/21)
- Kokkos needs to rever to ${PROJECT\_NAME}\_ENABLE\_CXX11 not Trilinos\_ENABLE\_CXX11 [\#17](https://github.com/kokkos/kokkos/issues/17)
- Kokkos Base Makefile adds AVX to KNC Build [\#16](https://github.com/kokkos/kokkos/issues/16)
- MS Visual Studio 2013 Build Errors [\#9](https://github.com/kokkos/kokkos/issues/9)
- subview\(X, ALL\(\), j\) for 2-D LayoutRight View X: should it view a column? [\#5](https://github.com/kokkos/kokkos/issues/5)

## [End_C++98](https://github.com/kokkos/kokkos/tree/End_C++98) (2015-04-15)


\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*
