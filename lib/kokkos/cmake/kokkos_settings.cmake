########################## NOTES ###############################################
# This files goal is to take CMake options found in kokkos_options.cmake but 
# possibly set from elsewhere 
#   (see: trilinos/cmake/ProjectCOmpilerPostConfig.cmake) 
# using CMake idioms and map them onto the KOKKOS_SETTINGS variables that gets 
# passed to the kokkos makefile configuration:
#  make -f ${CMAKE_SOURCE_DIR}/core/src/Makefile ${KOKKOS_SETTINGS} build-makefile-cmake-kokkos
# that generates KokkosCore_config.h and kokkos_generated_settings.cmake
# To understand how to form KOKKOS_SETTINGS, see
#     <KOKKOS_PATH>/Makefile.kokkos

#-------------------------------------------------------------------------------
#------------------------------- GENERAL OPTIONS -------------------------------
#-------------------------------------------------------------------------------

# Ensure that KOKKOS_ARCH is in the ARCH_LIST
foreach(arch ${KOKKOS_ARCH})
  list(FIND KOKKOS_ARCH_LIST ${arch} indx)
  if (indx EQUAL -1)
    message(FATAL_ERROR "${arch} is not an accepted value for KOKKOS_ARCH."
      "  Please pick from these choices: ${KOKKOS_INTERNAL_ARCH_DOCSTR}")
  endif ()
endforeach()

# KOKKOS_SETTINGS uses KOKKOS_ARCH
string(REPLACE ";" "," KOKKOS_ARCH "${KOKKOS_ARCH}")
set(KOKKOS_ARCH ${KOKKOS_ARCH})

# From Makefile.kokkos: Options: yes,no
if(${KOKKOS_ENABLE_DEBUG})
  set(KOKKOS_DEBUG yes)
else()
  set(KOKKOS_DEBUG no)
endif()

#------------------------------- KOKKOS_DEVICES --------------------------------
# Can have multiple devices 
set(KOKKOS_DEVICESl)
foreach(devopt ${KOKKOS_DEVICES_LIST})
  string(TOUPPER ${devopt} devoptuc)
  if (${KOKKOS_ENABLE_${devoptuc}}) 
    list(APPEND KOKKOS_DEVICESl ${devopt})
  endif ()
endforeach()
# List needs to be comma-delmitted
string(REPLACE ";" "," KOKKOS_DEVICES "${KOKKOS_DEVICESl}")

#------------------------------- KOKKOS_OPTIONS --------------------------------
# From Makefile.kokkos: Options: aggressive_vectorization,disable_profiling
#compiler_warnings, aggressive_vectorization, disable_profiling, disable_dualview_modify_check, enable_profile_load_print

set(KOKKOS_OPTIONSl)
if(${KOKKOS_ENABLE_COMPILER_WARNINGS})
      list(APPEND KOKKOS_OPTIONSl compiler_warnings)
endif()
if(${KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION})
      list(APPEND KOKKOS_OPTIONSl aggressive_vectorization)
endif()
if(NOT ${KOKKOS_ENABLE_PROFILING})
      list(APPEND KOKKOS_OPTIONSl disable_vectorization)
endif()
if(NOT ${KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK})
      list(APPEND KOKKOS_OPTIONSl disable_dualview_modify_check)
endif()
if(${KOKKOS_ENABLE_PROFILING_LOAD_PRINT})
      list(APPEND KOKKOS_OPTIONSl enable_profile_load_print)
endif()
# List needs to be comma-delimitted
string(REPLACE ";" "," KOKKOS_OPTIONS "${KOKKOS_OPTIONSl}")


#------------------------------- KOKKOS_USE_TPLS -------------------------------
# Construct the Makefile options
set(KOKKOS_USE_TPLSl)
foreach(tplopt ${KOKKOS_USE_TPLS_LIST})
  if (${KOKKOS_ENABLE_${tplopt}}) 
    list(APPEND KOKKOS_USE_TPLSl ${KOKKOS_INTERNAL_${tplopt}})
  endif ()
endforeach()
# List needs to be comma-delimitted
string(REPLACE ";" "," KOKKOS_USE_TPLS "${KOKKOS_USE_TPLSl}")


#------------------------------- KOKKOS_CUDA_OPTIONS ---------------------------
# Construct the Makefile options
set(KOKKOS_CUDA_OPTIONS)
foreach(cudaopt ${KOKKOS_CUDA_OPTIONS_LIST})
  if (${KOKKOS_ENABLE_CUDA_${cudaopt}})
    list(APPEND KOKKOS_CUDA_OPTIONSl ${KOKKOS_INTERNAL_${cudaopt}})
  endif ()
endforeach()
# List needs to be comma-delmitted
string(REPLACE ";" "," KOKKOS_CUDA_OPTIONS "${KOKKOS_CUDA_OPTIONSl}")

#------------------------------- PATH VARIABLES --------------------------------
#  Want makefile to use same executables specified which means modifying
#  the path so the $(shell ...) commands in the makefile see the right exec
#  Also, the Makefile's use FOO_PATH naming scheme for -I/-L construction
#TODO:  Makefile.kokkos allows this to be overwritten? ROCM_HCC_PATH

set(KOKKOS_INTERNAL_PATHS)
set(addpathl)
foreach(kvar "CUDA;QTHREADS;${KOKKOS_USE_TPLS_LIST}")
  if(${KOKKOS_ENABLE_${kvar}})
    if(DEFINED KOKKOS_${kvar}_DIR)
      set(KOKKOS_INTERNAL_PATHS "${KOKKOS_INTERNAL_PATHS} ${kvar}_PATH=${KOKKOS_${kvar}_DIR}")
      if(IS_DIRECTORY ${KOKKOS_${kvar}_DIR}/bin)
        list(APPEND addpathl ${KOKKOS_${kvar}_DIR}/bin)
      endif()
    endif()
  endif()
endforeach()
# Path env is : delimitted
string(REPLACE ";" ":" KOKKOS_INTERNAL_ADDTOPATH "${addpathl}")


######################### SET KOKKOS_SETTINGS ##################################
# Set the KOKKOS_SETTINGS String -- this is the primary communication with the
# makefile configuration.  See Makefile.kokkos

set(KOKKOS_SETTINGS KOKKOS_SRC_PATH=${KOKKOS_SRC_PATH})
set(KOKKOS_SETTINGS ${KOKKOS_SETTINGS} KOKKOS_PATH=${KOKKOS_PATH})
set(KOKKOS_SETTINGS ${KOKKOS_SETTINGS} KOKKOS_INSTALL_PATH=${CMAKE_INSTALL_PREFIX})

# Form of KOKKOS_foo=$KOKKOS_foo
foreach(kvar ARCH;DEVICES;DEBUG;OPTIONS;CUDA_OPTIONS;USE_TPLS)
  set(KOKKOS_VAR KOKKOS_${kvar})
  if(DEFINED KOKKOS_${kvar})
    if (NOT "${${KOKKOS_VAR}}" STREQUAL "")
      set(KOKKOS_SETTINGS ${KOKKOS_SETTINGS} ${KOKKOS_VAR}=${${KOKKOS_VAR}})
    endif()
  endif()
endforeach()

# Form of VAR=VAL
#TODO:  Makefile supports MPICH_CXX, OMPI_CXX as well
foreach(ovar CXX;CXXFLAGS;LDFLAGS)
  if(DEFINED ${ovar})
    if (NOT "${${ovar}}" STREQUAL "")
      set(KOKKOS_SETTINGS ${KOKKOS_SETTINGS} ${ovar}=${${ovar}})
    endif()
  endif()
endforeach()

# Finally, do the paths
if (NOT "${KOKKOS_INTERNAL_PATHS}" STREQUAL "")
  set(KOKKOS_SETTINGS ${KOKKOS_SETTINGS} ${KOKKOS_INTERNAL_PATHS})
endif()
if (NOT "${KOKKOS_INTERNAL_ADDTOPATH}" STREQUAL "")
  set(KOKKOS_SETTINGS ${KOKKOS_SETTINGS} PATH=${KOKKOS_INTERNAL_ADDTOPATH}:\${PATH})
endif()

# Final form that gets passed to make
set(KOKKOS_SETTINGS env ${KOKKOS_SETTINGS})


############################ PRINT CONFIGURE STATUS ############################

if(KOKKOS_CMAKE_VERBOSE)
  message(STATUS "")
  message(STATUS "****************** Kokkos Settings ******************")
  message(STATUS "Execution Spaces")

  if(KOKKOS_ENABLE_CUDA)
    message(STATUS "  Device Parallel: Cuda")
  else()
    message(STATUS "  Device Parallel: None")
  endif()

  if(KOKKOS_ENABLE_OPENMP)
    message(STATUS "    Host Parallel: OpenMP")
  elseif(KOKKOS_ENABLE_PTHREAD)
    message(STATUS "    Host Parallel: Pthread")
  elseif(KOKKOS_ENABLE_QTHREADS)
    message(STATUS "    Host Parallel: Qthreads")
  else()
    message(STATUS "    Host Parallel: None")
  endif()

  if(KOKKOS_ENABLE_SERIAL)
    message(STATUS "      Host Serial: Serial")
  else()
    message(STATUS "      Host Serial: None")
  endif()

  message(STATUS "")
  message(STATUS "Architectures:")
  message(STATUS "    ${KOKKOS_ARCH}")

  message(STATUS "")
  message(STATUS "Enabled options")

  if(KOKKOS_SEPARATE_LIBS)
    message(STATUS "  KOKKOS_SEPARATE_LIBS")
  endif()

  if(KOKKOS_ENABLE_HWLOC)
    message(STATUS "  KOKKOS_ENABLE_HWLOC")
  endif()

  if(KOKKOS_ENABLE_MEMKIND)
    message(STATUS "  KOKKOS_ENABLE_MEMKIND")
  endif()

  if(KOKKOS_ENABLE_DEBUG)
    message(STATUS "  KOKKOS_ENABLE_DEBUG")
  endif()

  if(KOKKOS_ENABLE_PROFILING)
    message(STATUS "  KOKKOS_ENABLE_PROFILING")
  endif()

  if(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION)
    message(STATUS "  KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION")
  endif()

  if(KOKKOS_ENABLE_CUDA)
    if(KOKKOS_ENABLE_CUDA_LDG_INTRINSIC)
      message(STATUS "  KOKKOS_ENABLE_CUDA_LDG_INTRINSIC")
    endif()

    if(KOKKOS_ENABLE_CUDA_UVM)
      message(STATUS "  KOKKOS_ENABLE_CUDA_UVM")
    endif()

    if(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
      message(STATUS "  KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE")
    endif()

    if(KOKKOS_ENABLE_CUDA_LAMBDA)
      message(STATUS "  KOKKOS_ENABLE_CUDA_LAMBDA")
    endif()

    if(KOKKOS_CUDA_DIR)
      message(STATUS "  KOKKOS_CUDA_DIR: ${KOKKOS_CUDA_DIR}")
    endif()
  endif()

  if(KOKKOS_QTHREADS_DIR)
    message(STATUS "  KOKKOS_QTHREADS_DIR: ${KOKKOS_QTHREADS_DIR}")
  endif()

  if(KOKKOS_HWLOC_DIR)
    message(STATUS "  KOKKOS_HWLOC_DIR: ${KOKKOS_HWLOC_DIR}")
  endif()

  if(KOKKOS_MEMKIND_DIR)
    message(STATUS "  KOKKOS_MEMKIND_DIR: ${KOKKOS_MEMKIND_DIR}")
  endif()

  message(STATUS "")
  message(STATUS "Final kokkos settings variable:")
  message(STATUS "  ${KOKKOS_SETTINGS}")

  message(STATUS "*****************************************************")
  message(STATUS "")
endif()
