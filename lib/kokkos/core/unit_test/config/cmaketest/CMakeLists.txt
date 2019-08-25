cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(Kokkos CXX)

enable_testing()

# Initialization
get_filename_component(KOKKOS_TESTDIR ${CMAKE_SOURCE_DIR}/../.. REALPATH)
get_filename_component(KOKKOS_SRCDIR ${CMAKE_SOURCE_DIR}/../../../.. REALPATH)
set(KOKKOS_SRC_PATH ${KOKKOS_SRCDIR})
set(KOKKOS_PATH ${KOKKOS_SRC_PATH})

set(CXX ${KOKKOS_TESTDIR}/config/cxx)

# Defined in core/src/Makefile -- this should be consistent
set(KOKKOS_MAKEFILE Makefile.kokkos)
set(KOKKOS_CMAKEFILE kokkos_generated_settings.cmake)

# Defined in Makefile.kokkos -- this should be consistent
set(KOKKOS_INTERNAL_CONFIG_TMP KokkosCore_config.tmp)
set(KOKKOS_CONFIG_HEADER KokkosCore_config.h)

set(KOKKOS_CMAKE_VERBOSE False)
include(${KOKKOS_SRCDIR}/cmake/kokkos_options.cmake)
foreach(KOKKOS_DEV ${KOKKOS_DEVICES_LIST})
# Do some initialization: Want to turn everything off for testing
  string(TOUPPER ${KOKKOS_DEV} KOKKOS_DEVUC)
  set(KOKKOS_ENABLE_${KOKKOS_DEVUC} OFF)
endforeach()


#TEST set(KOKKOS_HOST_ARCH_LIST ARMv80)
#TEST set(KOKKOS_DEVICES_LIST Cuda)
#set(KOKKOS_HOST_ARCH_LIST AMDAVX)
#set(KOKKOS_DEVICES_LIST Cuda)

foreach(KOKKOS_HOST_ARCH ${KOKKOS_HOST_ARCH_LIST})
  foreach(KOKKOS_DEV ${KOKKOS_DEVICES_LIST})
    string(TOUPPER ${KOKKOS_DEV} KOKKOS_DEVUC)
    set(KOKKOS_ENABLE_${KOKKOS_DEVUC} On)

    set(KOKKOS_CMAKE_VERBOSE True)
    include(${KOKKOS_SRCDIR}/cmake/kokkos_options.cmake)
    set(KOKKOS_SETTINGS ${KOKKOS_SETTINGS} ROCM_HCC_PATH=${KOKKOS_TESTDIR}/config)

    #message(STATUS "${KOKKOS_SETTINGS} make -f ${KOKKOS_SRCDIR}/core/src/Makefile build-makefile-cmake-kokkos")
    execute_process(
          COMMAND ${KOKKOS_SETTINGS} make -f ${KOKKOS_SRCDIR}/core/src/Makefile build-makefile-cmake-kokkos
          WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
          OUTPUT_FILE ${CMAKE_BINARY_DIR}/core_src_make.out
          RESULT_VARIABLE res
        )
    #message(STATUS "RESULT ${res}")

    file(REMOVE ${KOKKOS_INTERNAL_CONFIG_TMP} ${KOKKOS_MAKEFILE})
    set(PREFIX "${KOKKOS_HOST_ARCH}_${KOKKOS_DEV}_")
    set(NEWCMAKE ${PREFIX}${KOKKOS_CMAKEFILE})
    set(NEWCONFH ${PREFIX}${KOKKOS_CONFIG_HEADER})
    file(RENAME ${KOKKOS_CMAKEFILE} ${NEWCMAKE})
    file(RENAME ${KOKKOS_CONFIG_HEADER} ${NEWCONFH})
      
    add_test(NAME ${NEWCMAKE}-test 
        COMMAND ${KOKKOS_TESTDIR}/testmake.sh ${NEWCMAKE}  ${KOKKOS_HOST_ARCH} ${KOKKOS_DEV}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )
    set_tests_properties(${NEWCMAKE}-test 
        PROPERTIES PASS_REGULAR_EXPRESSION Passed
        TIMEOUT 15
      )
   add_test(NAME ${NEWCONFH}-test 
        COMMAND ${KOKKOS_TESTDIR}/diffconfig.sh ${NEWCONFH} 
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )
   set_tests_properties(${NEWCONFH}-test 
        PROPERTIES PASS_REGULAR_EXPRESSION Passed
        TIMEOUT 15
      )
    set(KOKKOS_ENABLE_${KOKKOS_DEVUC} Off)

  endforeach()
endforeach()
