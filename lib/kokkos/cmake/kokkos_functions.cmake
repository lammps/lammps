################################### FUNCTIONS ##################################
# List of functions
#   set_kokkos_cxx_compiler
#   set_kokkos_cxx_standard
#   set_kokkos_srcs

#-------------------------------------------------------------------------------
# function(set_kokkos_cxx_compiler)
# Sets the following compiler variables that are analogous to the CMAKE_*
# versions.  We add the ability to detect NVCC (really nvcc_wrapper).
#   KOKKOS_CXX_COMPILER
#   KOKKOS_CXX_COMPILER_ID
#   KOKKOS_CXX_COMPILER_VERSION
#
# Inputs:
#   KOKKOS_ENABLE_CUDA
#   CMAKE_CXX_COMPILER
#   CMAKE_CXX_COMPILER_ID
#   CMAKE_CXX_COMPILER_VERSION
#
# Also verifies the compiler version meets the minimum required by Kokkos.
function(set_kokkos_cxx_compiler)
  # Since CMake doesn't recognize the nvcc compiler until 3.8, we use our own
  # version of the CMake variables and detect nvcc ourselves.  Initially set to
  # the CMake variable values.
  set(INTERNAL_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  set(INTERNAL_CXX_COMPILER_ID ${CMAKE_CXX_COMPILER_ID})
  set(INTERNAL_CXX_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})

  # Check if the compiler is nvcc (which really means nvcc_wrapper).
  execute_process(COMMAND ${INTERNAL_CXX_COMPILER} --version
                  COMMAND grep nvcc
                  COMMAND wc -l
                  OUTPUT_VARIABLE INTERNAL_HAVE_COMPILER_NVCC
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  string(REGEX REPLACE "^ +" ""
         INTERNAL_HAVE_COMPILER_NVCC ${INTERNAL_HAVE_COMPILER_NVCC})

  if(INTERNAL_HAVE_COMPILER_NVCC)
    # Set the compiler id to nvcc.  We use the value used by CMake 3.8.
    set(INTERNAL_CXX_COMPILER_ID NVIDIA)

    # Set nvcc's compiler version.
    execute_process(COMMAND ${INTERNAL_CXX_COMPILER} --version
                    COMMAND grep release
                    OUTPUT_VARIABLE INTERNAL_CXX_COMPILER_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+$"
           INTERNAL_CXX_COMPILER_VERSION ${INTERNAL_CXX_COMPILER_VERSION})
  endif()

  # Enforce the minimum compilers supported by Kokkos.
  set(KOKKOS_MESSAGE_TEXT "Compiler not supported by Kokkos.  Required compiler versions:")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    Clang      3.5.2 or higher")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    GCC        4.8.4 or higher")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    Intel     15.0.2 or higher")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    NVCC      7.0.28 or higher")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    PGI         17.1 or higher\n")

  if(INTERNAL_CXX_COMPILER_ID STREQUAL Clang)
    if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 3.5.2)
      message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
    endif()
  elseif(INTERNAL_CXX_COMPILER_ID STREQUAL GNU)
    if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 4.8.4)
      message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
    endif()
  elseif(INTERNAL_CXX_COMPILER_ID STREQUAL Intel)
    if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 15.0.2)
      message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
    endif()
  elseif(INTERNAL_CXX_COMPILER_ID STREQUAL NVIDIA)
    if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 7.0.28)
      message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
    endif()
  elseif(INTERNAL_CXX_COMPILER_ID STREQUAL PGI)
    if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 17.1)
      message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
    endif()
  endif()

  # Enforce that extensions are turned off for nvcc_wrapper.
  if(INTERNAL_CXX_COMPILER_ID STREQUAL NVIDIA)
    if(DEFINED CMAKE_CXX_EXTENSIONS AND CMAKE_CXX_EXTENSIONS STREQUAL ON)
      message(FATAL_ERROR "NVCC doesn't support C++ extensions.  Set CMAKE_CXX_EXTENSIONS to OFF in your CMakeLists.txt.")
    endif()
  endif()

  if(KOKKOS_ENABLE_CUDA)
    # Enforce that the compiler can compile CUDA code.
    if(INTERNAL_CXX_COMPILER_ID STREQUAL Clang)
      if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 4.0.0)
        message(FATAL_ERROR "Compiling CUDA code directly with Clang requires version 4.0.0 or higher.")
      endif()
    elseif(NOT INTERNAL_CXX_COMPILER_ID STREQUAL NVIDIA)
      message(FATAL_ERROR "Invalid compiler for CUDA.  The compiler must be nvcc_wrapper or Clang.")
    endif()
  endif()

  set(KOKKOS_CXX_COMPILER ${INTERNAL_CXX_COMPILER} PARENT_SCOPE)
  set(KOKKOS_CXX_COMPILER_ID ${INTERNAL_CXX_COMPILER_ID} PARENT_SCOPE)
  set(KOKKOS_CXX_COMPILER_VERSION ${INTERNAL_CXX_COMPILER_VERSION} PARENT_SCOPE)
endfunction()

#-------------------------------------------------------------------------------
# function(set_kokkos_cxx_standard)
#  Transitively enforces that the appropriate CXX standard compile flags (C++11
#  or above) are added to targets that use the Kokkos library.  Compile features
#  are used if possible.  Otherwise, the appropriate flags are added to
#  KOKKOS_CXX_FLAGS.  Values set by the user to CMAKE_CXX_STANDARD and
#  CMAKE_CXX_EXTENSIONS are honored.
#
# Outputs:
#   KOKKOS_CXX11_FEATURES
#   KOKKOS_CXX_FLAGS
#
# Inputs:
#  KOKKOS_CXX_COMPILER
#  KOKKOS_CXX_COMPILER_ID
#  KOKKOS_CXX_COMPILER_VERSION
#
function(set_kokkos_cxx_standard)
  # The following table lists the versions of CMake that supports CXX_STANDARD
  # and the CXX compile features for different compilers.  The versions are
  # based on CMake documentation, looking at CMake code, and verifying by
  # testing with specific CMake versions.
  #
  #   COMPILER                      CXX_STANDARD     Compile Features
  #   ---------------------------------------------------------------
  #   Clang                             3.1                3.1
  #   GNU                               3.1                3.2
  #   AppleClang                        3.2                3.2
  #   Intel                             3.6                3.6
  #   Cray                              No                 No
  #   PGI                               No                 No
  #   XL                                No                 No
  #
  # For compiling CUDA code using nvcc_wrapper, we will use the host compiler's
  # flags for turning on C++11.  Since for compiler ID and versioning purposes
  # CMake recognizes the host compiler when calling nvcc_wrapper, this just
  # works.  Both NVCC and nvcc_wrapper only recognize '-std=c++11' which means
  # that we can only use host compilers for CUDA builds that use those flags.
  # It also means that extensions (gnu++11) can't be turned on for CUDA builds.

  # Check if we can use compile features.
  if(NOT KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
    if(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
      if(NOT CMAKE_VERSION VERSION_LESS 3.1)
        set(INTERNAL_USE_COMPILE_FEATURES ON)
      endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL AppleClang OR CMAKE_CXX_COMPILER_ID STREQUAL GNU)
      if(NOT CMAKE_VERSION VERSION_LESS 3.2)
        set(INTERNAL_USE_COMPILE_FEATURES ON)
      endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      if(NOT CMAKE_VERSION VERSION_LESS 3.6)
        set(INTERNAL_USE_COMPILE_FEATURES ON)
      endif()
    endif()
  endif()

  if(INTERNAL_USE_COMPILE_FEATURES)
    # Use the compile features aspect of CMake to transitively cause C++ flags
    # to populate to user code.

    # I'm using a hack by requiring features that I know force the lowest version
    # of the compilers we want to support.  Clang 3.3 and later support all of
    # the C++11 standard.  With CMake 3.8 and higher, we could switch to using
    # cxx_std_11.
    set(KOKKOS_CXX11_FEATURES
        cxx_nonstatic_member_init # Forces GCC 4.7 or later and Intel 14.0 or later.
        PARENT_SCOPE
       )
  else()
    # CXX compile features are not yet implemented for this combination of
    # compiler and version of CMake.

    if(CMAKE_CXX_COMPILER_ID STREQUAL AppleClang)
      # Versions of CMAKE before 3.2 don't support CXX_STANDARD or C++ compile
      # features for the AppleClang compiler.  Set compiler flags transitively
      # here such that they trickle down to a call to target_compile_options().

      # The following two blocks of code were copied from
      # /Modules/Compiler/AppleClang-CXX.cmake from CMake 3.7.2 and then
      # modified.
      if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0)
        set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
        set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++11")
      endif()

      if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1)
        set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "-std=c++14")
        set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++14")
      elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.1)
        # AppleClang 5.0 knows this flag, but does not set a __cplusplus macro
        # greater than 201103L.
        set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "-std=c++1y")
        set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++1y")
      endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      # Versions of CMAKE before 3.6 don't support CXX_STANDARD or C++ compile
      # features for the Intel compiler.  Set compiler flags transitively here
      # such that they trickle down to a call to target_compile_options().

      # The following three blocks of code were copied from
      # /Modules/Compiler/Intel-CXX.cmake from CMake 3.7.2 and then modified.
      if("x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
        set(_std -Qstd)
        set(_ext c++)
      else()
        set(_std -std)
        set(_ext gnu++)
      endif()

      if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.2)
        set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "${_std}=c++14")
        # TODO: There is no gnu++14 value supported; figure out what to do.
        set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "${_std}=c++14")
      elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.0)
        set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "${_std}=c++1y")
        # TODO: There is no gnu++14 value supported; figure out what to do.
        set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "${_std}=c++1y")
      endif()

      if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0)
        set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "${_std}=c++11")
        set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "${_std}=${_ext}11")
      elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.1)
        set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "${_std}=c++0x")
        set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "${_std}=${_ext}0x")
      endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL Cray)
      # CMAKE doesn't support CXX_STANDARD or C++ compile features for the Cray
      # compiler.  Set compiler options transitively here such that they trickle
      # down to a call to target_compile_options().
      set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "-hstd=c++11")
      set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "-hstd=c++11")
      set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "-hstd=c++11")
      set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "-hstd=c++11")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL PGI)
      # CMAKE doesn't support CXX_STANDARD or C++ compile features for the PGI
      # compiler.  Set compiler options transitively here such that they trickle
      # down to a call to target_compile_options().
      set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "--c++11")
      set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "--c++11")
      set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "--c++11")
      set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "--c++11")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL XL)
      # CMAKE doesn't support CXX_STANDARD or C++ compile features for the XL
      # compiler.  Set compiler options transitively here such that they trickle
      # down to a call to target_compile_options().
      set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
      set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "-std=c++11")
      set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "-std=c++11")
      set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "-std=c++11")
    else()
      # Assume GNU.  CMAKE_CXX_STANDARD is handled correctly by CMake 3.1 and
      # above for this compiler.  If the user explicitly requests a C++
      # standard, CMake takes care of it.  If not, transitively require C++11.
      if(NOT CMAKE_CXX_STANDARD)
        set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION ${CMAKE_CXX11_STANDARD_COMPILE_OPTION})
        set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION ${CMAKE_CXX11_EXTENSION_COMPILE_OPTION})
      endif()
    endif()

    # Set the C++ standard info for Kokkos respecting user set values for
    # CMAKE_CXX_STANDARD and CMAKE_CXX_EXTENSIONS.
    # Only use cxx extension if explicitly requested
    if(CMAKE_CXX_STANDARD EQUAL 14)
      if(DEFINED CMAKE_CXX_EXTENSIONS AND CMAKE_CXX_EXTENSIONS STREQUAL ON)
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX14_EXTENSION_COMPILE_OPTION})
      else()
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX14_STANDARD_COMPILE_OPTION})
      endif()
    elseif(CMAKE_CXX_STANDARD EQUAL 11)
      if(DEFINED CMAKE_CXX_EXTENSIONS AND CMAKE_CXX_EXTENSIONS STREQUAL ON)
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX11_EXTENSION_COMPILE_OPTION})
      else()
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX11_STANDARD_COMPILE_OPTION})
      endif()
    else()
      # The user didn't explicitly request a standard, transitively require
      # C++11 respecting CMAKE_CXX_EXTENSIONS.
      if(DEFINED CMAKE_CXX_EXTENSIONS AND CMAKE_CXX_EXTENSIONS STREQUAL ON)
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX11_EXTENSION_COMPILE_OPTION})
      else()
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX11_STANDARD_COMPILE_OPTION})
      endif()
    endif()

    set(KOKKOS_CXX_FLAGS ${INTERNAL_CXX_FLAGS} PARENT_SCOPE)
  endif()
endfunction()


#-------------------------------------------------------------------------------
# function(set_kokkos_sources)
# Takes a list of sources for kokkos (e.g., KOKKOS_SRC from Makefile.kokkos and
# put it into kokkos_generated_settings.cmake) and sorts the files into the subpackages or
# separate_libraries.  This is core and containers (algorithms is pure header
# files).
#
# Inputs:
#   KOKKOS_SRC
# 
# Outputs:
#   KOKKOS_CORE_SRCS
#   KOKKOS_CONTAINERS_SRCS
#
function(set_kokkos_srcs)
  set(opts ) # no-value args
  set(oneValArgs )
  set(multValArgs KOKKOS_SRC) # e.g., lists
  cmake_parse_arguments(IN "${opts}" "${oneValArgs}" "${multValArgs}" ${ARGN})

  foreach(sfile ${IN_KOKKOS_SRC})
     string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" stripfile "${sfile}")
     string(REPLACE "/" ";" striplist "${stripfile}")
     list(GET striplist 0 firstdir)
     if(${firstdir} STREQUAL "core")
       list(APPEND KOKKOS_CORE_SRCS ${sfile})
     else()
       list(APPEND KOKKOS_CONTAINERS_SRCS ${sfile})
     endif()
  endforeach()
  set(KOKKOS_CORE_SRCS ${KOKKOS_CORE_SRCS} PARENT_SCOPE)
  set(KOKKOS_CONTAINERS_SRCS ${KOKKOS_CONTAINERS_SRCS} PARENT_SCOPE)
  return()
endfunction()

# Setting a default value if it is not already set
macro(set_kokkos_default_default VARIABLE DEFAULT)
  IF( "${KOKKOS_INTERNAL_ENABLE_${VARIABLE}_DEFAULT}" STREQUAL "" )
    IF( "${KOKKOS_ENABLE_${VARIABLE}}" STREQUAL "" )
      set(KOKKOS_INTERNAL_ENABLE_${VARIABLE}_DEFAULT ${DEFAULT})
  #    MESSAGE(WARNING "Set: KOKKOS_INTERNAL_ENABLE_${VARIABLE}_DEFAULT to ${KOKKOS_INTERNAL_ENABLE_${VARIABLE}_DEFAULT}")
    ELSE()
      set(KOKKOS_INTERNAL_ENABLE_${VARIABLE}_DEFAULT ${KOKKOS_ENABLE_${VARIABLE}})
   #   MESSAGE(WARNING "Set: KOKKOS_INTERNAL_ENABLE_${VARIABLE}_DEFAULT to ${KOKKOS_INTERNAL_ENABLE_${VARIABLE}_DEFAULT}")
    ENDIF()
  ENDIF()
  UNSET(KOKKOS_ENABLE_${VARIABLE} CACHE)
endmacro()
