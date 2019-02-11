########################## NOTES ###############################################
#  List the options for configuring kokkos using CMake method of doing it.
#  These options then get mapped onto KOKKOS_SETTINGS environment variable by
#  kokkos_settings.cmake.  It is separate to allow other packages to override
#  these variables (e.g., TriBITS).

########################## AVAILABLE OPTIONS ###################################
# Use lists for documentation, verification, and programming convenience

# All CMake options of the type KOKKOS_ENABLE_*
set(KOKKOS_INTERNAL_ENABLE_OPTIONS_LIST)
list(APPEND KOKKOS_INTERNAL_ENABLE_OPTIONS_LIST
     Serial
     OpenMP
     Pthread
     Qthread
     Cuda
     ROCm
     HWLOC
     MEMKIND
     LIBRT
     Cuda_Lambda
     Cuda_Relocatable_Device_Code
     Cuda_UVM
     Cuda_LDG_Intrinsic
     Debug
     Debug_DualView_Modify_Check
     Debug_Bounds_Check
     Compiler_Warnings
     Profiling
     Profiling_Load_Print
     Aggressive_Vectorization
     Deprecated_Code
     Explicit_Instantiation
     )

#-------------------------------------------------------------------------------
#------------------------------- Recognize CamelCase Options ---------------------------
#-------------------------------------------------------------------------------

foreach(opt ${KOKKOS_INTERNAL_ENABLE_OPTIONS_LIST})
  string(TOUPPER ${opt} OPT )
  IF(DEFINED Kokkos_ENABLE_${opt})
    IF(DEFINED KOKKOS_ENABLE_${OPT})
      IF(NOT ("${KOKKOS_ENABLE_${OPT}}" STREQUAL "${Kokkos_ENABLE_${opt}}"))
        IF(DEFINED KOKKOS_ENABLE_${OPT}_INTERNAL)
          MESSAGE(WARNING  "Defined both Kokkos_ENABLE_${opt}=[${Kokkos_ENABLE_${opt}}] and KOKKOS_ENABLE_${OPT}=[${KOKKOS_ENABLE_${OPT}}] and they differ! Could be caused by old CMakeCache Variable. Run CMake again and warning should disappear. If not you are truly setting both variables.")
          IF(NOT ("${Kokkos_ENABLE_${opt}}" STREQUAL "${KOKKOS_ENABLE_${OPT}_INTERNAL}"))
            UNSET(KOKKOS_ENABLE_${OPT} CACHE)
            SET(KOKKOS_ENABLE_${OPT} ${Kokkos_ENABLE_${opt}})
            MESSAGE(WARNING "SET BOTH VARIABLES KOKKOS_ENABLE_${OPT}: ${KOKKOS_ENABLE_${OPT}}")
          ELSE()
            SET(Kokkos_ENABLE_${opt} ${KOKKOS_ENABLE_${OPT}})
          ENDIF()
        ELSE()
          MESSAGE(FATAL_ERROR "Defined both Kokkos_ENABLE_${opt}=[${Kokkos_ENABLE_${opt}}] and KOKKOS_ENABLE_${OPT}=[${KOKKOS_ENABLE_${OPT}}] and they differ!")
        ENDIF()
      ENDIF()
    ELSE()
      SET(KOKKOS_INTERNAL_ENABLE_${OPT}_DEFAULT ${Kokkos_ENABLE_${opt}})
    ENDIF()
  ENDIF()
endforeach()

IF(DEFINED Kokkos_ARCH)
  MESSAGE(FATAL_ERROR "Defined Kokkos_ARCH, use KOKKOS_ARCH instead!")
ENDIF()
IF(DEFINED Kokkos_Arch)
  MESSAGE(FATAL_ERROR "Defined Kokkos_Arch, use KOKKOS_ARCH instead!")
ENDIF()
  
#-------------------------------------------------------------------------------
# List of possible host architectures.
#-------------------------------------------------------------------------------
set(KOKKOS_ARCH_LIST)
list(APPEND KOKKOS_ARCH_LIST
     None            # No architecture optimization
     AMDAVX          # (HOST) AMD chip
     ARMv80          # (HOST) ARMv8.0 Compatible CPU
     ARMv81          # (HOST) ARMv8.1 Compatible CPU
     ARMv8-ThunderX  # (HOST) ARMv8 Cavium ThunderX CPU
     ARMv8-TX2       # (HOST) ARMv8 Cavium ThunderX2 CPU
     WSM             # (HOST) Intel Westmere CPU
     SNB             # (HOST) Intel Sandy/Ivy Bridge CPUs
     HSW             # (HOST) Intel Haswell CPUs
     BDW             # (HOST) Intel Broadwell Xeon E-class CPUs
     SKX             # (HOST) Intel Sky Lake Xeon E-class HPC CPUs (AVX512)
     KNC             # (HOST) Intel Knights Corner Xeon Phi
     KNL             # (HOST) Intel Knights Landing Xeon Phi
     BGQ             # (HOST) IBM Blue Gene Q
     Power7          # (HOST) IBM POWER7 CPUs
     Power8          # (HOST) IBM POWER8 CPUs
     Power9          # (HOST) IBM POWER9 CPUs
     Kepler          # (GPU) NVIDIA Kepler default (generation CC 3.5)
     Kepler30        # (GPU) NVIDIA Kepler generation CC 3.0
     Kepler32        # (GPU) NVIDIA Kepler generation CC 3.2
     Kepler35        # (GPU) NVIDIA Kepler generation CC 3.5
     Kepler37        # (GPU) NVIDIA Kepler generation CC 3.7
     Maxwell         # (GPU) NVIDIA Maxwell default (generation CC 5.0)
     Maxwell50       # (GPU) NVIDIA Maxwell generation CC 5.0
     Maxwell52       # (GPU) NVIDIA Maxwell generation CC 5.2
     Maxwell53       # (GPU) NVIDIA Maxwell generation CC 5.3
     Pascal60        # (GPU) NVIDIA Pascal generation CC 6.0
     Pascal61        # (GPU) NVIDIA Pascal generation CC 6.1
     Volta70         # (GPU) NVIDIA Volta generation CC 7.0
     Volta72         # (GPU) NVIDIA Volta generation CC 7.2
     Turing75         # (GPU) NVIDIA Turing generation CC 7.5
    )

# List of possible device architectures.
# The case and spelling here needs to match Makefile.kokkos
set(KOKKOS_DEVICES_LIST)
# Options: Cuda,ROCm,OpenMP,Pthread,Qthreads,Serial
list(APPEND KOKKOS_DEVICES_LIST
    Cuda          # NVIDIA GPU -- see below
    OpenMP        # OpenMP
    Pthread       # pthread
    Qthreads      # qthreads
    Serial        # serial
    ROCm          # Relocatable device code
    )

# List of possible TPLs for Kokkos
# From Makefile.kokkos: Options: hwloc,librt,experimental_memkind
set(KOKKOS_USE_TPLS_LIST)
if(APPLE)
list(APPEND KOKKOS_USE_TPLS_LIST
    HWLOC          # hwloc
    MEMKIND        # experimental_memkind
    )
else()
list(APPEND KOKKOS_USE_TPLS_LIST
    HWLOC          # hwloc
    LIBRT          # librt
    MEMKIND        # experimental_memkind
    )
endif()
# Map of cmake variables to Makefile variables
set(KOKKOS_INTERNAL_HWLOC hwloc)
set(KOKKOS_INTERNAL_LIBRT librt)
set(KOKKOS_INTERNAL_MEMKIND experimental_memkind)

# List of possible Advanced options
set(KOKKOS_OPTIONS_LIST)
list(APPEND KOKKOS_OPTIONS_LIST
       AGGRESSIVE_VECTORIZATION    
       DISABLE_PROFILING          
       DISABLE_DUALVIEW_MODIFY_CHECK
       ENABLE_PROFILE_LOAD_PRINT   
    )
# Map of cmake variables to Makefile variables
set(KOKKOS_INTERNAL_LDG_INTRINSIC use_ldg)
set(KOKKOS_INTERNAL_UVM librt)
set(KOKKOS_INTERNAL_RELOCATABLE_DEVICE_CODE rdc)


#-------------------------------------------------------------------------------
# List of possible Options for CUDA
#-------------------------------------------------------------------------------
# From Makefile.kokkos: Options: use_ldg,force_uvm,rdc
set(KOKKOS_CUDA_OPTIONS_LIST)
list(APPEND KOKKOS_CUDA_OPTIONS_LIST
    LDG_INTRINSIC              # use_ldg
    UVM                        # force_uvm
    RELOCATABLE_DEVICE_CODE    # rdc
    LAMBDA                     # enable_lambda
    )
    
# Map of cmake variables to Makefile variables
set(KOKKOS_INTERNAL_LDG_INTRINSIC use_ldg)
set(KOKKOS_INTERNAL_UVM force_uvm)
set(KOKKOS_INTERNAL_RELOCATABLE_DEVICE_CODE rdc)
set(KOKKOS_INTERNAL_LAMBDA enable_lambda)


#-------------------------------------------------------------------------------
#------------------------------- Create doc strings ----------------------------
#-------------------------------------------------------------------------------

set(tmpr "\n       ")
string(REPLACE ";" ${tmpr} KOKKOS_INTERNAL_ARCH_DOCSTR "${KOKKOS_ARCH_LIST}")
set(KOKKOS_INTERNAL_ARCH_DOCSTR "${tmpr}${KOKKOS_INTERNAL_ARCH_DOCSTR}")
# This would be useful, but we use Foo_ENABLE mechanisms
#string(REPLACE ";" ${tmpr} KOKKOS_INTERNAL_DEVICES_DOCSTR "${KOKKOS_DEVICES_LIST}")
#string(REPLACE ";" ${tmpr} KOKKOS_INTERNAL_USE_TPLS_DOCSTR "${KOKKOS_USE_TPLS_LIST}")
#string(REPLACE ";" ${tmpr} KOKKOS_INTERNAL_CUDA_OPTIONS_DOCSTR "${KOKKOS_CUDA_OPTIONS_LIST}")

#-------------------------------------------------------------------------------
#------------------------------- GENERAL OPTIONS -------------------------------
#-------------------------------------------------------------------------------

# Setting this variable to a value other than "None" can improve host
# performance by turning on architecture specific code.
# NOT SET is used to determine if the option is passed in.  It is reset to
# default "None" down below.
set(KOKKOS_ARCH "NOT_SET" CACHE STRING 
      "Optimize for specific host architecture. Options are: ${KOKKOS_INTERNAL_ARCH_DOCSTR}")

# Whether to build separate libraries or now
set(KOKKOS_SEPARATE_LIBS OFF CACHE BOOL "OFF = kokkos.  ON = kokkoscore, kokkoscontainers, and kokkosalgorithms.")

# Qthreads options.
set(KOKKOS_QTHREADS_DIR "" CACHE PATH "Location of Qthreads library.")


#-------------------------------------------------------------------------------
#------------------------------- KOKKOS_DEVICES --------------------------------
#-------------------------------------------------------------------------------
# Figure out default settings
IF(Trilinos_ENABLE_Kokkos)             
  set_kokkos_default_default(SERIAL ON)
  set_kokkos_default_default(PTHREAD OFF)
  IF(TPL_ENABLE_QTHREAD)
    set_kokkos_default_default(QTHREADS ${TPL_ENABLE_QTHREAD})
  ELSE()
    set_kokkos_default_default(QTHREADS OFF)
  ENDIF()
  IF(Trilinos_ENABLE_OpenMP)
    set_kokkos_default_default(OPENMP ${Trilinos_ENABLE_OpenMP})
  ELSE()
    set_kokkos_default_default(OPENMP OFF)
  ENDIF()
  IF(TPL_ENABLE_CUDA)
    set_kokkos_default_default(CUDA ${TPL_ENABLE_CUDA})
  ELSE()
    set_kokkos_default_default(CUDA OFF)
  ENDIF()
  set_kokkos_default_default(ROCM OFF)
ELSE()
  set_kokkos_default_default(SERIAL ON)
  set_kokkos_default_default(OPENMP OFF)
  set_kokkos_default_default(PTHREAD OFF)
  set_kokkos_default_default(QTHREAD OFF)
  set_kokkos_default_default(CUDA OFF)
  set_kokkos_default_default(ROCM OFF)
ENDIF()

# Set which Kokkos backend to use.
# These are the actual options that define the settings.
set(KOKKOS_ENABLE_SERIAL ${KOKKOS_INTERNAL_ENABLE_SERIAL_DEFAULT} CACHE BOOL "Whether to enable the Kokkos::Serial device.  This device executes \"parallel\" kernels sequentially on a single CPU thread.  It is enabled by default.  If you disable this device, please enable at least one other CPU device, such as Kokkos::OpenMP or Kokkos::Threads.")
set(KOKKOS_ENABLE_OPENMP ${KOKKOS_INTERNAL_ENABLE_OPENMP_DEFAULT} CACHE BOOL "Enable OpenMP support in Kokkos." FORCE)
set(KOKKOS_ENABLE_PTHREAD ${KOKKOS_INTERNAL_ENABLE_PTHREAD_DEFAULT} CACHE BOOL "Enable Pthread support in Kokkos.")
set(KOKKOS_ENABLE_QTHREADS ${KOKKOS_INTERNAL_ENABLE_QTHREADS_DEFAULT} CACHE BOOL "Enable Qthreads support in Kokkos.")
set(KOKKOS_ENABLE_CUDA ${KOKKOS_INTERNAL_ENABLE_CUDA_DEFAULT} CACHE BOOL "Enable CUDA support in Kokkos.")
set(KOKKOS_ENABLE_ROCM ${KOKKOS_INTERNAL_ENABLE_ROCM_DEFAULT} CACHE BOOL "Enable ROCm support in Kokkos.")



#-------------------------------------------------------------------------------
#------------------------------- KOKKOS DEBUG and PROFILING --------------------
#-------------------------------------------------------------------------------

# Debug related options enable compiler warnings

set_kokkos_default_default(DEBUG OFF)
set(KOKKOS_ENABLE_DEBUG ${KOKKOS_INTERNAL_ENABLE_DEBUG_DEFAULT} CACHE BOOL "Enable Kokkos Debug.")

# From Makefile.kokkos: Advanced Options: 
#compiler_warnings, aggressive_vectorization, disable_profiling, disable_dualview_modify_check, enable_profile_load_print
set_kokkos_default_default(COMPILER_WARNINGS OFF)
set(KOKKOS_ENABLE_COMPILER_WARNINGS ${KOKKOS_INTERNAL_ENABLE_COMPILER_WARNINGS_DEFAULT} CACHE BOOL "Enable compiler warnings.")

set_kokkos_default_default(DEBUG_DUALVIEW_MODIFY_CHECK OFF)
set(KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK ${KOKKOS_INTERNAL_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK_DEFAULT} CACHE BOOL "Enable dualview modify check.")

# Enable aggressive vectorization.
set_kokkos_default_default(AGGRESSIVE_VECTORIZATION OFF)
set(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION ${KOKKOS_INTERNAL_ENABLE_AGGRESSIVE_VECTORIZATION_DEFAULT} CACHE BOOL "Enable aggressive vectorization.")

# Enable profiling.
set_kokkos_default_default(PROFILING ON)
set(KOKKOS_ENABLE_PROFILING ${KOKKOS_INTERNAL_ENABLE_PROFILING_DEFAULT} CACHE BOOL "Enable profiling.")

set_kokkos_default_default(PROFILING_LOAD_PRINT OFF)
set(KOKKOS_ENABLE_PROFILING_LOAD_PRINT ${KOKKOS_INTERNAL_ENABLE_PROFILING_LOAD_PRINT_DEFAULT} CACHE BOOL "Enable profile load print.")

set_kokkos_default_default(DEPRECATED_CODE ON)
set(KOKKOS_ENABLE_DEPRECATED_CODE ${KOKKOS_INTERNAL_ENABLE_DEPRECATED_CODE_DEFAULT} CACHE BOOL "Enable deprecated code.")

set_kokkos_default_default(EXPLICIT_INSTANTIATION OFF)
set(KOKKOS_ENABLE_EXPLICIT_INSTANTIATION ${KOKKOS_INTERNAL_ENABLE_EXPLICIT_INSTANTIATION_DEFAULT} CACHE BOOL "Enable explicit template instantiation.")

#-------------------------------------------------------------------------------
#------------------------------- KOKKOS_USE_TPLS -------------------------------
#-------------------------------------------------------------------------------
# Enable hwloc library.
# Figure out default:
IF(Trilinos_ENABLE_Kokkos AND TPL_ENABLE_HWLOC)
  set_kokkos_default_default(HWLOC ON)
ELSE()
  set_kokkos_default_default(HWLOC OFF)
ENDIF()
set(KOKKOS_ENABLE_HWLOC ${KOKKOS_INTERNAL_ENABLE_HWLOC_DEFAULT} CACHE BOOL "Enable hwloc for better process placement.")
set(KOKKOS_HWLOC_DIR "" CACHE PATH "Location of hwloc library. (kokkos tpl)")

# Enable memkind library.
set_kokkos_default_default(MEMKIND OFF)
set(KOKKOS_ENABLE_MEMKIND ${KOKKOS_INTERNAL_ENABLE_MEMKIND_DEFAULT} CACHE BOOL "Enable memkind. (kokkos tpl)")
set(KOKKOS_MEMKIND_DIR "" CACHE PATH "Location of memkind library. (kokkos tpl)")

# Enable rt library.
IF(Trilinos_ENABLE_Kokkos)
  IF(DEFINED TPL_ENABLE_LIBRT)
    set_kokkos_default_default(LIBRT ${TPL_ENABLE_LIBRT})
  ELSE()
    set_kokkos_default_default(LIBRT OFF)
  ENDIF()
ELSE()
  set_kokkos_default_default(LIBRT ON)
ENDIF()
set(KOKKOS_ENABLE_LIBRT ${KOKKOS_INTERNAL_ENABLE_LIBRT_DEFAULT} CACHE BOOL "Enable librt for more precise timer.  (kokkos tpl)")


#-------------------------------------------------------------------------------
#------------------------------- KOKKOS_CUDA_OPTIONS ---------------------------
#-------------------------------------------------------------------------------

# CUDA options.
# Set Defaults
set_kokkos_default_default(CUDA_LDG_INTRINSIC_DEFAULT OFF)
set_kokkos_default_default(CUDA_UVM_DEFAULT OFF)
set_kokkos_default_default(CUDA_RELOCATABLE_DEVICE_CODE OFF)
IF(Trilinos_ENABLE_Kokkos)
  IF(KOKKOS_ENABLE_CUDA)
    find_package(CUDA)
  ENDIF()
  IF (DEFINED CUDA_VERSION)
    IF (CUDA_VERSION VERSION_GREATER "7.0")
      set_kokkos_default_default(CUDA_LAMBDA ON)
    ELSE()
      set_kokkos_default_default(CUDA_LAMBDA OFF)
    ENDIF()
  ENDIF()
ELSE()
  set_kokkos_default_default(CUDA_LAMBDA OFF)
ENDIF()

# Set actual options
set(KOKKOS_CUDA_DIR "" CACHE PATH "Location of CUDA library.  Defaults to where nvcc installed.")
set(KOKKOS_ENABLE_CUDA_LDG_INTRINSIC ${KOKKOS_INTERNAL_ENABLE_CUDA_LDG_INTRINSIC_DEFAULT} CACHE BOOL "Enable CUDA LDG. (cuda option)") 
set(KOKKOS_ENABLE_CUDA_UVM ${KOKKOS_INTERNAL_ENABLE_CUDA_UVM_DEFAULT} CACHE BOOL "Enable CUDA unified virtual memory.")
set(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ${KOKKOS_INTERNAL_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE_DEFAULT} CACHE BOOL "Enable relocatable device code for CUDA. (cuda option)")
set(KOKKOS_ENABLE_CUDA_LAMBDA ${KOKKOS_INTERNAL_ENABLE_CUDA_LAMBDA_DEFAULT} CACHE BOOL "Enable lambdas for CUDA. (cuda option)")


#-------------------------------------------------------------------------------
#----------------------- HOST ARCH AND LEGACY TRIBITS --------------------------
#-------------------------------------------------------------------------------

# This defines the previous legacy TriBITS builds. 
set(KOKKOS_LEGACY_TRIBITS False)
IF ("${KOKKOS_ARCH}" STREQUAL "NOT_SET")
  set(KOKKOS_ARCH "None")
  IF(KOKKOS_HAS_TRILINOS)
    set(KOKKOS_LEGACY_TRIBITS True)
  ENDIF()
ENDIF()
IF (KOKKOS_HAS_TRILINOS)
  IF (KOKKOS_LEGACY_TRIBITS)
    message(STATUS "Using the legacy tribits build because KOKKOS_ARCH not set")
  ELSE()
    message(STATUS "NOT using the legacy tribits build because KOKKOS_ARCH *is* set")
  ENDIF()
ENDIF()

#-------------------------------------------------------------------------------
#----------------------- Set CamelCase Options if they are not yet set ---------
#-------------------------------------------------------------------------------

foreach(opt ${KOKKOS_INTERNAL_ENABLE_OPTIONS_LIST})
  string(TOUPPER ${opt} OPT )
  UNSET(KOKKOS_ENABLE_${OPT}_INTERNAL CACHE)
  SET(KOKKOS_ENABLE_${OPT}_INTERNAL ${KOKKOS_ENABLE_${OPT}} CACHE BOOL INTERNAL)
  IF(DEFINED KOKKOS_ENABLE_${OPT})
    UNSET(Kokkos_ENABLE_${opt} CACHE)
    SET(Kokkos_ENABLE_${opt} ${KOKKOS_ENABLE_${OPT}} CACHE BOOL "CamelCase Compatibility setting for KOKKOS_ENABLE_${OPT}")
  ENDIF()
endforeach()

