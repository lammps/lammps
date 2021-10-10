message(STATUS "Downloading and building OpenCL loader library")
set(OPENCL_LOADER_URL "${LAMMPS_THIRDPARTY_URL}/opencl-loader-2021.09.18.tar.gz" CACHE STRING "URL for OpenCL loader tarball")
set(OPENCL_LOADER_MD5 "3b3882627964bd02e5c3b02065daac3c" CACHE STRING "MD5 checksum of OpenCL loader tarball")
mark_as_advanced(OPENCL_LOADER_URL)
mark_as_advanced(OPENCL_LOADER_MD5)

include(ExternalProject)
ExternalProject_Add(opencl_loader
                    URL             ${OPENCL_LOADER_URL}
                    URL_MD5         ${OPENCL_LOADER_MD5}
                    SOURCE_DIR      "${CMAKE_BINARY_DIR}/opencl_loader-src"
                    BINARY_DIR      "${CMAKE_BINARY_DIR}/opencl_loader-build"
                    CMAKE_ARGS      ${CMAKE_REQUEST_PIC} ${CMAKE_EXTRA_OPENCL_LOADER_OPTS}
                                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                                    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                                    -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
                                    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
                    BUILD_BYPRODUCTS <BINARY_DIR>/libOpenCL${CMAKE_STATIC_LIBRARY_SUFFIX}
                    LOG_DOWNLOAD ON
                    LOG_CONFIGURE ON
                    LOG_BUILD ON
                    INSTALL_COMMAND ""
                    TEST_COMMAND    "")

ExternalProject_Get_Property(opencl_loader SOURCE_DIR)
set(OPENCL_LOADER_INCLUDE_DIR ${SOURCE_DIR}/inc)

# workaround for CMake 3.10 on ubuntu 18.04
file(MAKE_DIRECTORY ${OPENCL_LOADER_INCLUDE_DIR})

ExternalProject_Get_Property(opencl_loader BINARY_DIR)
set(OPENCL_LOADER_LIBRARY_PATH "${BINARY_DIR}/libOpenCL${CMAKE_STATIC_LIBRARY_SUFFIX}")

find_package(Threads QUIET)
if(NOT WIN32)
  set(OPENCL_LOADER_DEP_LIBS "Threads::Threads;${CMAKE_DL_LIBS}")
else()
  set(OPENCL_LOADER_DEP_LIBS "cfgmgr32;runtimeobject")
endif()

add_library(OpenCL::OpenCL UNKNOWN IMPORTED)
add_dependencies(OpenCL::OpenCL opencl_loader)

set_target_properties(OpenCL::OpenCL PROPERTIES
  IMPORTED_LOCATION ${OPENCL_LOADER_LIBRARY_PATH}
  INTERFACE_INCLUDE_DIRECTORIES ${OPENCL_LOADER_INCLUDE_DIR}
  INTERFACE_LINK_LIBRARIES "${OPENCL_LOADER_DEP_LIBS}")


