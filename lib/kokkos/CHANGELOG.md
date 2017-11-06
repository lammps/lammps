# Change Log
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
