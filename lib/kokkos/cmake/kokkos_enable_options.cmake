########################## NOTES ###############################################
#  List the options for configuring kokkos using CMake method of doing it.
#  These options then get mapped onto KOKKOS_SETTINGS environment variable by
#  kokkos_settings.cmake.  It is separate to allow other packages to override
#  these variables (e.g., TriBITS).

########################## AVAILABLE OPTIONS ###################################
# Use lists for documentation, verification, and programming convenience


FUNCTION(KOKKOS_ENABLE_OPTION SUFFIX DEFAULT DOCSTRING)
  KOKKOS_OPTION(ENABLE_${SUFFIX} ${DEFAULT} BOOL ${DOCSTRING})
  STRING(TOUPPER ${SUFFIX} UC_NAME)
  IF (KOKKOS_ENABLE_${UC_NAME} AND NOT "Kokkos_ENABLE_${UC_NAME}" IN_LIST Kokkos_OPTIONS_NOT_TO_EXPORT)
    LIST(APPEND KOKKOS_ENABLED_OPTIONS ${UC_NAME})
    #I hate that CMake makes me do this
    SET(KOKKOS_ENABLED_OPTIONS ${KOKKOS_ENABLED_OPTIONS} PARENT_SCOPE)
  ENDIF()
  SET(KOKKOS_ENABLE_${UC_NAME} ${KOKKOS_ENABLE_${UC_NAME}} PARENT_SCOPE)
ENDFUNCTION()

# Certain defaults will depend on knowing the enabled devices
KOKKOS_CFG_DEPENDS(OPTIONS DEVICES)
KOKKOS_CFG_DEPENDS(OPTIONS COMPILER_ID)

# Put a check in just in case people are using this option
KOKKOS_DEPRECATED_LIST(OPTIONS ENABLE)

KOKKOS_ENABLE_OPTION(CUDA_RELOCATABLE_DEVICE_CODE  OFF "Whether to enable relocatable device code (RDC) for CUDA")
KOKKOS_ENABLE_OPTION(CUDA_UVM             OFF "Whether to use unified memory (UM) for CUDA by default")
KOKKOS_ENABLE_OPTION(CUDA_LDG_INTRINSIC   OFF "Whether to use CUDA LDG intrinsics")
# In contrast to other CUDA-dependent, options CUDA_LAMBDA is ON by default.
# That is problematic when CUDA is not enabled because this not only yields a
# bogus warning, but also exports the Kokkos_ENABLE_CUDA_LAMBDA variable and
# sets it to ON. This if-clause is a crutch that delays the refactoring of the
# way we declare all options until after we get rid of TriBITS.
IF (Trilinos_ENABLE_Kokkos AND TPL_ENABLE_CUDA)
   SET(CUDA_LAMBDA_DEFAULT ON)
ELSEIF (KOKKOS_ENABLE_CUDA)
   SET(CUDA_LAMBDA_DEFAULT ON)
ELSE()
   SET(CUDA_LAMBDA_DEFAULT OFF)
ENDIF()
KOKKOS_ENABLE_OPTION(CUDA_LAMBDA ${CUDA_LAMBDA_DEFAULT} "Whether to allow lambda expressions on the device with NVCC **DEPRECATED**")

# May be used to disable our use of CudaMallocAsync.  It had caused issues in
# the past when UCX was used as MPI communication layer.  We expect it is
# resolved but we keep the option around a bit longer to be safe.
KOKKOS_ENABLE_OPTION(IMPL_CUDA_MALLOC_ASYNC ON  "Whether to enable CudaMallocAsync (requires CUDA Toolkit 11.2)")
KOKKOS_ENABLE_OPTION(IMPL_NVHPC_AS_DEVICE_COMPILER OFF "Whether to allow nvc++ as Cuda device compiler")
KOKKOS_ENABLE_OPTION(DEPRECATED_CODE_3    OFF "Whether code deprecated in major release 3 is available" )
KOKKOS_ENABLE_OPTION(DEPRECATED_CODE_4    ON "Whether code deprecated in major release 4 is available" )
KOKKOS_ENABLE_OPTION(DEPRECATION_WARNINGS ON "Whether to emit deprecation warnings" )
KOKKOS_ENABLE_OPTION(HIP_RELOCATABLE_DEVICE_CODE  OFF "Whether to enable relocatable device code (RDC) for HIP")
KOKKOS_ENABLE_OPTION(TESTS         OFF  "Whether to build the unit tests")
KOKKOS_ENABLE_OPTION(BENCHMARKS    OFF  "Whether to build the benchmarks")
KOKKOS_ENABLE_OPTION(EXAMPLES      OFF  "Whether to build the examples")
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" UPPERCASE_CMAKE_BUILD_TYPE)
IF(UPPERCASE_CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  KOKKOS_ENABLE_OPTION(DEBUG                ON "Whether to activate extra debug features - may increase compile times")
  KOKKOS_ENABLE_OPTION(DEBUG_DUALVIEW_MODIFY_CHECK ON "Debug check on dual views")
ELSE()
  KOKKOS_ENABLE_OPTION(DEBUG                OFF "Whether to activate extra debug features - may increase compile times")
  KOKKOS_ENABLE_OPTION(DEBUG_DUALVIEW_MODIFY_CHECK OFF "Debug check on dual views")
ENDIF()
UNSET(_UPPERCASE_CMAKE_BUILD_TYPE)
KOKKOS_ENABLE_OPTION(LARGE_MEM_TESTS      OFF "Whether to perform extra large memory tests")
KOKKOS_ENABLE_OPTION(DEBUG_BOUNDS_CHECK   OFF "Whether to use bounds checking - will increase runtime")
KOKKOS_ENABLE_OPTION(COMPILER_WARNINGS    OFF "Whether to print all compiler warnings")
KOKKOS_ENABLE_OPTION(TUNING               OFF "Whether to create bindings for tuning tools")
KOKKOS_ENABLE_OPTION(AGGRESSIVE_VECTORIZATION OFF "Whether to aggressively vectorize loops")
KOKKOS_ENABLE_OPTION(COMPILE_AS_CMAKE_LANGUAGE OFF "Whether to use native cmake language support")
KOKKOS_ENABLE_OPTION(HIP_MULTIPLE_KERNEL_INSTANTIATIONS OFF "Whether multiple kernels are instantiated at compile time - improve performance but increase compile time")

# This option will go away eventually, but allows fallback to old implementation when needed.
KOKKOS_ENABLE_OPTION(DESUL_ATOMICS_EXTERNAL OFF "Whether to use an external desul installation")

KOKKOS_ENABLE_OPTION(IMPL_MDSPAN OFF "Whether to enable experimental mdspan support")
KOKKOS_ENABLE_OPTION(MDSPAN_EXTERNAL OFF BOOL "Whether to use an external version of mdspan")
KOKKOS_ENABLE_OPTION(IMPL_SKIP_COMPILER_MDSPAN ON BOOL "Whether to use an internal version of mdspan even if the compiler supports mdspan")
mark_as_advanced(Kokkos_ENABLE_IMPL_MDSPAN)
mark_as_advanced(Kokkos_ENABLE_MDSPAN_EXTERNAL)
mark_as_advanced(Kokkos_ENABLE_IMPL_SKIP_COMPILER_MDSPAN)

IF (Trilinos_ENABLE_Kokkos)
  SET(COMPLEX_ALIGN_DEFAULT OFF)
ELSE()
  SET(COMPLEX_ALIGN_DEFAULT ON)
ENDIF()
KOKKOS_ENABLE_OPTION(COMPLEX_ALIGN ${COMPLEX_ALIGN_DEFAULT}  "Whether to align Kokkos::complex to 2*alignof(RealType)")

IF (KOKKOS_ENABLE_TESTS)
  SET(HEADER_SELF_CONTAINMENT_TESTS_DEFAULT ON)
ELSE()
  SET(HEADER_SELF_CONTAINMENT_TESTS_DEFAULT OFF)
ENDIF()
KOKKOS_ENABLE_OPTION(HEADER_SELF_CONTAINMENT_TESTS ${HEADER_SELF_CONTAINMENT_TESTS_DEFAULT} "Enable header self-containment unit tests")
IF (NOT KOKKOS_ENABLE_TESTS AND KOKKOS_ENABLE_HEADER_SELF_CONTAINMENT_TESTS)
  MESSAGE(WARNING "Kokkos_ENABLE_HEADER_SELF_CONTAINMENT_TESTS is ON but Kokkos_ENABLE_TESTS is OFF. Option will be ignored.")
ENDIF()

IF (KOKKOS_ENABLE_CUDA AND (KOKKOS_CXX_COMPILER_ID STREQUAL Clang))
  SET(CUDA_CONSTEXPR_DEFAULT ON)
ELSE()
  SET(CUDA_CONSTEXPR_DEFAULT OFF)
ENDIF()
KOKKOS_ENABLE_OPTION(CUDA_CONSTEXPR ${CUDA_CONSTEXPR_DEFAULT} "Whether to activate experimental relaxed constexpr functions")

IF (KOKKOS_ENABLE_HPX)
  SET(HPX_ASYNC_DISPATCH_DEFAULT ON)
ELSE()
  SET(HPX_ASYNC_DISPATCH_DEFAULT OFF)
ENDIF()
KOKKOS_ENABLE_OPTION(IMPL_HPX_ASYNC_DISPATCH ${HPX_ASYNC_DISPATCH_DEFAULT} "Whether HPX supports asynchronous dispatch")

Kokkos_ENABLE_OPTION(UNSUPPORTED_ARCHS OFF "Whether to allow architectures in backends Kokkos doesn't optimize for")

FUNCTION(check_device_specific_options)
  CMAKE_PARSE_ARGUMENTS(SOME "" "DEVICE" "OPTIONS" ${ARGN})
  IF(NOT KOKKOS_ENABLE_${SOME_DEVICE})
    FOREACH(OPTION ${SOME_OPTIONS})
      IF(NOT DEFINED CACHE{Kokkos_ENABLE_${OPTION}} OR NOT DEFINED CACHE{Kokkos_ENABLE_${SOME_DEVICE}})
        MESSAGE(FATAL_ERROR "Internal logic error: option '${OPTION}' or device '${SOME_DEVICE}' not recognized.")
      ENDIF()
      IF(KOKKOS_ENABLE_${OPTION})
        MESSAGE(WARNING "Kokkos_ENABLE_${OPTION} is ON but ${SOME_DEVICE} backend is not enabled. Option will be ignored.")
        UNSET(KOKKOS_ENABLE_${OPTION} PARENT_SCOPE)
      ENDIF()
    ENDFOREACH()
  ENDIF()
ENDFUNCTION()

CHECK_DEVICE_SPECIFIC_OPTIONS(DEVICE CUDA OPTIONS CUDA_UVM CUDA_RELOCATABLE_DEVICE_CODE CUDA_LAMBDA CUDA_CONSTEXPR CUDA_LDG_INTRINSIC)
CHECK_DEVICE_SPECIFIC_OPTIONS(DEVICE HIP OPTIONS HIP_RELOCATABLE_DEVICE_CODE)
CHECK_DEVICE_SPECIFIC_OPTIONS(DEVICE HPX OPTIONS IMPL_HPX_ASYNC_DISPATCH)

# Needed due to change from deprecated name to new header define name
IF (KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION)
  SET(KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION ON)
ENDIF()

# Force consistency of KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE
# and CMAKE_CUDA_SEPARABLE_COMPILATION when we are compiling
# using the CMake CUDA language support.
# Either one being on will turn the other one on.
IF (KOKKOS_COMPILE_LANGUAGE STREQUAL CUDA)
  IF (KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
    IF (NOT CMAKE_CUDA_SEPARABLE_COMPILATION)
      MESSAGE(STATUS "Setting CMAKE_CUDA_SEPARABLE_COMPILATION=ON since Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE is true. When compiling Kokkos with CMake language CUDA, please use CMAKE_CUDA_SEPARABLE_COMPILATION to control RDC support")
      SET(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
    ENDIF()
  ELSE()
    IF (CMAKE_CUDA_SEPARABLE_COMPILATION)
      SET(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON)
    ENDIF()
  ENDIF()
ENDIF()

# This is known to occur with Clang 9. We would need to use nvcc as the linker
# http://lists.llvm.org/pipermail/cfe-dev/2018-June/058296.html
# TODO: Through great effort we can use a different linker by hacking
# CMAKE_CXX_LINK_EXECUTABLE in a future release
IF (KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE AND KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
  MESSAGE(FATAL_ERROR "Relocatable device code is currently not supported with Clang - must use nvcc_wrapper or turn off RDC")
ENDIF()

IF (KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE AND BUILD_SHARED_LIBS)
  MESSAGE(FATAL_ERROR "Relocatable device code requires static libraries.")
ENDIF()

IF(Kokkos_ENABLE_CUDA_LDG_INTRINSIC)
  IF(KOKKOS_ENABLE_DEPRECATED_CODE_4)
    MESSAGE(DEPRECATION "Setting Kokkos_ENABLE_CUDA_LDG_INTRINSIC is deprecated. LDG intrinsics are always enabled.")
  ELSE()
    MESSAGE(FATAL_ERROR "Kokkos_ENABLE_CUDA_LDG_INTRINSIC has been removed. LDG intrinsics are always enabled.")
  ENDIF()
ENDIF()
IF(Kokkos_ENABLE_CUDA AND NOT Kokkos_ENABLE_CUDA_LAMBDA)
  IF(KOKKOS_ENABLE_DEPRECATED_CODE_4)
    MESSAGE(DEPRECATION "Setting Kokkos_ENABLE_CUDA_LAMBDA is deprecated. Lambda expressions in device code are always enabled. Forcing -DKokkos_ENABLE_CUDA_LAMBDA=ON")
    set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Kokkos turned Cuda lambda support ON!" FORCE)
    set(KOKKOS_ENABLE_CUDA_LAMBDA ON)
  ELSE()
    MESSAGE(FATAL_ERROR "Kokkos_ENABLE_CUDA_LAMBDA has been removed. Lambda expressions in device code always enabled.")
  ENDIF()
ENDIF()


IF(DEFINED Kokkos_ENABLE_IMPL_DESUL_ATOMICS)
  MESSAGE(WARNING "Kokkos_ENABLE_IMPL_DESUL_ATOMICS option has been removed. Desul atomics cannot be disabled.")
ENDIF()
