KOKKOS_CFG_DEPENDS(COMPILER_ID NONE)

SET(KOKKOS_CXX_COMPILER ${CMAKE_CXX_COMPILER})
SET(KOKKOS_CXX_COMPILER_ID ${CMAKE_CXX_COMPILER_ID})
SET(KOKKOS_CXX_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})

MACRO(kokkos_internal_have_compiler_nvcc)
  # Check if the compiler is nvcc (which really means nvcc_wrapper).
  EXECUTE_PROCESS(COMMAND ${ARGN} --version
                  OUTPUT_VARIABLE INTERNAL_COMPILER_VERSION
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  STRING(REPLACE "\n" " - " INTERNAL_COMPILER_VERSION_ONE_LINE ${INTERNAL_COMPILER_VERSION} )
  STRING(FIND ${INTERNAL_COMPILER_VERSION_ONE_LINE} "nvcc" INTERNAL_COMPILER_VERSION_CONTAINS_NVCC)
  STRING(REGEX REPLACE "^ +" "" INTERNAL_HAVE_COMPILER_NVCC "${INTERNAL_HAVE_COMPILER_NVCC}")
  IF(${INTERNAL_COMPILER_VERSION_CONTAINS_NVCC} GREATER -1)
    SET(INTERNAL_HAVE_COMPILER_NVCC true)
  ELSE()
    SET(INTERNAL_HAVE_COMPILER_NVCC false)
  ENDIF()
ENDMACRO()

IF(Kokkos_ENABLE_CUDA)
  # find kokkos_launch_compiler
  FIND_PROGRAM(Kokkos_COMPILE_LAUNCHER
      NAMES           kokkos_launch_compiler
      HINTS           ${PROJECT_SOURCE_DIR}
      PATHS           ${PROJECT_SOURCE_DIR}
      PATH_SUFFIXES   bin)

  # check if compiler was set to nvcc_wrapper
  kokkos_internal_have_compiler_nvcc(${CMAKE_CXX_COMPILER})
  # if launcher was found and nvcc_wrapper was not specified as
  # compiler, set to use launcher. Will ensure CMAKE_CXX_COMPILER
  # is replaced by nvcc_wrapper
  IF(Kokkos_COMPILE_LAUNCHER AND NOT INTERNAL_HAVE_COMPILER_NVCC AND NOT KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
    # the first argument to launcher is always the C++ compiler defined by cmake
    # if the second argument matches the C++ compiler, it forwards the rest of the
    # args to nvcc_wrapper
    kokkos_internal_have_compiler_nvcc(
      ${Kokkos_COMPILE_LAUNCHER} ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER} -DKOKKOS_DEPENDENCE)
    SET(INTERNAL_USE_COMPILER_LAUNCHER true)
  ENDIF()
ENDIF()

IF(INTERNAL_HAVE_COMPILER_NVCC)
  # Save the host compiler id before overwriting it.
  SET(KOKKOS_CXX_HOST_COMPILER_ID ${KOKKOS_CXX_COMPILER_ID})

  # SET the compiler id to nvcc.  We use the value used by CMake 3.8.
  SET(KOKKOS_CXX_COMPILER_ID NVIDIA CACHE STRING INTERNAL FORCE)

  STRING(REGEX MATCH "V[0-9]+\\.[0-9]+\\.[0-9]+"
         TEMP_CXX_COMPILER_VERSION ${INTERNAL_COMPILER_VERSION_ONE_LINE})
  STRING(SUBSTRING ${TEMP_CXX_COMPILER_VERSION} 1 -1 TEMP_CXX_COMPILER_VERSION)
  SET(KOKKOS_CXX_COMPILER_VERSION ${TEMP_CXX_COMPILER_VERSION} CACHE STRING INTERNAL FORCE)
  MESSAGE(STATUS "Compiler Version: ${KOKKOS_CXX_COMPILER_VERSION}")
  IF(INTERNAL_USE_COMPILER_LAUNCHER)
    IF(Kokkos_LAUNCH_COMPILER_INFO)
        GET_FILENAME_COMPONENT(BASE_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
        # does not have STATUS intentionally
        MESSAGE("")
        MESSAGE("Kokkos_LAUNCH_COMPILER_INFO (${Kokkos_COMPILE_LAUNCHER}):")
        MESSAGE("  - Kokkos + CUDA backend requires the C++ files to be compiled as CUDA code.")
        MESSAGE("  - kokkos_launch_compiler permits CMAKE_CXX_COMPILER to be set to a traditional C++ compiler when Kokkos_ENABLE_CUDA=ON")
        MESSAGE("    by prefixing all the compile and link commands with the path to the script + CMAKE_CXX_COMPILER (${CMAKE_CXX_COMPILER}).")
        MESSAGE("  - If any of the compile or link commands have CMAKE_CXX_COMPILER as the first argument, it replaces CMAKE_CXX_COMPILER with nvcc_wrapper.")
        MESSAGE("  - If the compile or link command is not CMAKE_CXX_COMPILER, it just executes the command.")
        MESSAGE("  - If using ccache, set CMAKE_CXX_COMPILER to nvcc_wrapper explicitly.")
        MESSAGE("  - kokkos_compiler_launcher is available to downstream projects as well.")
        MESSAGE("    - If CMAKE_CXX_COMPILER=nvcc_wrapper, all legacy behavior will be preserved during 'find_package(Kokkos)'")
        MESSAGE("    - If CMAKE_CXX_COMPILER is not nvcc_wrapper, 'find_package(Kokkos)' will apply 'kokkos_compilation(GLOBAL)' unless separable compilation is enabled")
        MESSAGE("      - This can be disabled via '-DKokkos_LAUNCH_COMPILER=OFF'")
        MESSAGE("    - Use 'find_package(Kokkos COMPONENTS separable_compilation)' to enable separable compilation")
        MESSAGE("      - Separable compilation allows you to control the scope of where the compiler transformation behavior (${BASE_COMPILER_NAME} -> nvcc_wrapper) is applied")
        MESSAGE("      - The compiler transformation can be applied on a per-project, per-directory, per-target, and/or per-source-file basis")
        MESSAGE("        - 'kokkos_compilation(PROJECT)' will apply the compiler transformation to all targets in a project/subproject")
        MESSAGE("        - 'kokkos_compilation(TARGET <TARGET> [<TARGETS>...])' will apply the compiler transformation to the specified target(s)")
        MESSAGE("        - 'kokkos_compilation(SOURCE <SOURCE> [<SOURCES>...])' will apply the compiler transformation to the specified source file(s)")
        MESSAGE("        - 'kokkos_compilation(DIRECTORY <DIR> [<DIRS>...])' will apply the compiler transformation to the specified directories")
        MESSAGE("")
    ELSE()
        MESSAGE(STATUS "kokkos_launch_compiler (${Kokkos_COMPILE_LAUNCHER}) is enabled... Set Kokkos_LAUNCH_COMPILER_INFO=ON for more info.")
    ENDIF()
    kokkos_compilation(GLOBAL)
  ENDIF()
ENDIF()

IF(Kokkos_ENABLE_HIP)
  # get HIP version
  EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} --version
                  OUTPUT_VARIABLE INTERNAL_COMPILER_VERSION
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  STRING(REPLACE "\n" " - " INTERNAL_COMPILER_VERSION_ONE_LINE ${INTERNAL_COMPILER_VERSION} )
  SET(KOKKOS_CXX_COMPILER_ID HIP CACHE STRING INTERNAL FORCE)

  STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
         TEMP_CXX_COMPILER_VERSION ${INTERNAL_COMPILER_VERSION_ONE_LINE})
  SET(KOKKOS_CXX_COMPILER_VERSION ${TEMP_CXX_COMPILER_VERSION} CACHE STRING INTERNAL FORCE)
  MESSAGE(STATUS "Compiler Version: ${KOKKOS_CXX_COMPILER_VERSION}")
ENDIF()

IF(KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
  # The Cray compiler reports as Clang to most versions of CMake
  EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} --version
                  COMMAND grep Cray
                  COMMAND wc -l
                  OUTPUT_VARIABLE INTERNAL_HAVE_CRAY_COMPILER
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  IF (INTERNAL_HAVE_CRAY_COMPILER) #not actually Clang
    SET(KOKKOS_CLANG_IS_CRAY TRUE)
  ENDIF()
  # The clang based Intel compiler reports as Clang to most versions of CMake
  EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} --version
                  COMMAND grep icpx
                  COMMAND wc -l
                  OUTPUT_VARIABLE INTERNAL_HAVE_INTEL_COMPILER
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  IF (INTERNAL_HAVE_INTEL_COMPILER) #not actually Clang
    SET(KOKKOS_CLANG_IS_INTEL TRUE)
    SET(KOKKOS_CXX_COMPILER_ID IntelClang CACHE STRING INTERNAL FORCE)
  ENDIF()
ENDIF()

IF(KOKKOS_CXX_COMPILER_ID STREQUAL Cray OR KOKKOS_CLANG_IS_CRAY)
  # SET Cray's compiler version.
  EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} --version
                  OUTPUT_VARIABLE INTERNAL_CXX_COMPILER_VERSION
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
         TEMP_CXX_COMPILER_VERSION ${INTERNAL_CXX_COMPILER_VERSION})
  IF (KOKKOS_CLANG_IS_CRAY)
    SET(KOKKOS_CLANG_CRAY_COMPILER_VERSION ${TEMP_CXX_COMPILER_VERSION})
  ELSE()
    SET(KOKKOS_CXX_COMPILER_VERSION ${TEMP_CXX_COMPILER_VERSION} CACHE STRING INTERNAL FORCE)
  ENDIF()
ENDIF()

IF(KOKKOS_CXX_COMPILER_ID STREQUAL Fujitsu)
  # SET Fujitsus compiler version which is not detected by CMake
  EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} --version
                  OUTPUT_VARIABLE INTERNAL_CXX_COMPILER_VERSION
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
         TEMP_CXX_COMPILER_VERSION ${INTERNAL_CXX_COMPILER_VERSION})
  SET(KOKKOS_CXX_COMPILER_VERSION ${TEMP_CXX_COMPILER_VERSION} CACHE STRING INTERNAL FORCE)
ENDIF()

# Enforce the minimum compilers supported by Kokkos.
SET(KOKKOS_MESSAGE_TEXT "Compiler not supported by Kokkos.  Required compiler versions:")
SET(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    Clang      4.0.0 or higher")
SET(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    GCC        5.3.0 or higher")
SET(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    Intel     17.0.0 or higher")
SET(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    NVCC      9.2.88 or higher")
SET(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    HIPCC      3.8.0 or higher")
SET(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    PGI         17.4 or higher\n")

IF(KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
  IF(KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 4.0.0)
    MESSAGE(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
  ENDIF()
ELSEIF(KOKKOS_CXX_COMPILER_ID STREQUAL GNU)
  IF(KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 5.3.0)
    MESSAGE(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
  ENDIF()
ELSEIF(KOKKOS_CXX_COMPILER_ID STREQUAL Intel)
  IF(KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 17.0.0)
    MESSAGE(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
  ENDIF()
ELSEIF(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
  IF(KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 9.2.88)
    MESSAGE(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
  ENDIF()
  SET(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Kokkos turns off CXX extensions" FORCE)
ELSEIF(KOKKOS_CXX_COMPILER_ID STREQUAL HIP)
  IF(KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 3.8.0)
    MESSAGE(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
  ENDIF()
ELSEIF(KOKKOS_CXX_COMPILER_ID STREQUAL PGI)
  IF(KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 17.4)
    MESSAGE(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
  ENDIF()
ENDIF()

STRING(REPLACE "." ";" VERSION_LIST ${KOKKOS_CXX_COMPILER_VERSION})
LIST(GET VERSION_LIST 0 KOKKOS_COMPILER_VERSION_MAJOR)
LIST(GET VERSION_LIST 1 KOKKOS_COMPILER_VERSION_MINOR)
LIST(GET VERSION_LIST 2 KOKKOS_COMPILER_VERSION_PATCH)
