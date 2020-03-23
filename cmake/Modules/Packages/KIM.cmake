if(PKG_KIM)
  set(KIM-API_MIN_VERSION 2.1)
  find_package(CURL)
  if(CURL_FOUND)
    target_link_libraries(lammps PRIVATE CURL::libcurl) 
    add_definitions(-DLMP_KIM_CURL)
    set(LMP_DEBUG_CURL OFF CACHE STRING "Set libcurl verbose mode on/off. If on, it displays a lot of verbose information about its operations.")
    mark_as_advanced(LMP_DEBUG_CURL)
    if(LMP_DEBUG_CURL)
      add_definitions(-DLMP_DEBUG_CURL)
    endif()
    set(LMP_NO_SSL_CHECK OFF CACHE STRING "Tell libcurl to not verify the peer. If on, the connection succeeds regardless of the names in the certificate. Insecure - Use with caution!")
    mark_as_advanced(LMP_NO_SSL_CHECK)
    if(LMP_NO_SSL_CHECK)
      add_definitions(-DLMP_NO_SSL_CHECK)
    endif()
  endif()
  find_package(KIM-API QUIET)
  if(KIM-API_FOUND)
    if (KIM-API_VERSION VERSION_LESS ${KIM-API_MIN_VERSION})
      if ("${DOWNLOAD_KIM}" STREQUAL "")
        message(WARNING "Unsuitable KIM-API version \"${KIM-API_VERSION}\" found, but required is at least \"${KIM-API_MIN_VERSION}\".  Default behavior set to download and build our own.")
      endif()
      set(DOWNLOAD_KIM_DEFAULT ON)
    else()
      set(DOWNLOAD_KIM_DEFAULT OFF)
    endif()
  else()
    if ("${DOWNLOAD_KIM}" STREQUAL "")
      message(WARNING "KIM-API package not found.  Default behavior set to download and build our own")
    endif()
    set(DOWNLOAD_KIM_DEFAULT ON)
  endif()
  option(DOWNLOAD_KIM "Download KIM-API from OpenKIM instead of using an already installed one" ${DOWNLOAD_KIM_DEFAULT})
  if(DOWNLOAD_KIM)
    message(STATUS "KIM-API download requested - we will build our own")
    # Workaround for cross compilation with MinGW where ${CMAKE_INSTALL_LIBDIR}
    # is a full path, so we need to remove the prefix
    string(REPLACE ${CMAKE_INSTALL_PREFIX} "" _KIM_LIBDIR ${CMAKE_INSTALL_LIBDIR})
    include(ExternalProject)
    enable_language(C)
    enable_language(Fortran)
    ExternalProject_Add(kim_build
      URL https://s3.openkim.org/kim-api/kim-api-2.1.3.txz
      URL_MD5 6ee829a1bbba5f8b9874c88c4c4ebff8
      BINARY_DIR build
      CMAKE_ARGS ${CMAKE_REQUEST_PIC}
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                 -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                 -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                 -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
                 -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
                 BUILD_BYPRODUCTS <INSTALL_DIR>/${_KIM_LIBDIR}/libkim-api${CMAKE_SHARED_LIBRARY_SUFFIX} 
      )
    ExternalProject_get_property(kim_build INSTALL_DIR)
    set(KIM-API_INCLUDE_DIRS ${INSTALL_DIR}/include/kim-api)
    set(KIM-API_LDFLAGS ${INSTALL_DIR}/${_KIM_LIBDIR}/libkim-api${CMAKE_SHARED_LIBRARY_SUFFIX})
    add_dependencies(lammps kim_build)
  else()
    find_package(KIM-API ${KIM-API_MIN_VERSION} REQUIRED)
  endif()
  target_link_libraries(lammps PRIVATE "${KIM-API_LDFLAGS}")
  include_directories(${KIM-API_INCLUDE_DIRS})
endif()
