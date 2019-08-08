if(PKG_KIM)
  find_package(CURL)
  if(CURL_FOUND)
    include_directories(${CURL_INCLUDE_DIRS})
    list(APPEND LAMMPS_LINK_LIBS ${CURL_LIBRARIES})
    add_definitions(-DLMP_KIM_CURL)
  endif()
  find_package(KIM-API QUIET)
  if(KIM-API_FOUND)
    set(DOWNLOAD_KIM_DEFAULT OFF)
  else()
    if (NOT DOWNLOAD_KIM)
      message(WARNING "KIM-API package not found.  We will download and build our own")
    endif()
    set(DOWNLOAD_KIM_DEFAULT ON)
  endif()
  option(DOWNLOAD_KIM "Download KIM-API from OpenKIM instead of using an already installed one" ${DOWNLOAD_KIM_DEFAULT})
  if(DOWNLOAD_KIM)
    if(CMAKE_GENERATOR STREQUAL "Ninja")
      message(FATAL_ERROR "Cannot build downloaded KIM-API library with Ninja build tool")
    endif()
    message(STATUS "KIM-API download requested - we will build our own")
    include(CheckLanguage)
    include(ExternalProject)
    enable_language(C)
    check_language(Fortran)
    if(NOT CMAKE_Fortran_COMPILER)
      message(FATAL_ERROR "Compiling the KIM-API library requires a Fortran compiler")
    endif()
    ExternalProject_Add(kim_build
      URL https://s3.openkim.org/kim-api/kim-api-2.1.2.txz
      URL_MD5 6ac52e14ef52967fc7858220b208cba5
      BINARY_DIR build
      CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                 -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                 -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      )
    ExternalProject_get_property(kim_build INSTALL_DIR)
    set(KIM-API_INCLUDE_DIRS ${INSTALL_DIR}/include/kim-api)
    set(KIM-API_LDFLAGS ${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libkim-api${CMAKE_SHARED_LIBRARY_SUFFIX})
    list(APPEND LAMMPS_DEPS kim_build)
  else()
    find_package(KIM-API REQUIRED)
  endif()
  list(APPEND LAMMPS_LINK_LIBS "${KIM-API_LDFLAGS}")
  include_directories(${KIM-API_INCLUDE_DIRS})
endif()
