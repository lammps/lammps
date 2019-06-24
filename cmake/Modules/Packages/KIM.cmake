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
    set(DOWNLOAD_KIM_DEFAULT ON)
  endif()
  option(DOWNLOAD_KIM "Download KIM-API from OpenKIM instead of using an already installed one" ${DOWNLOAD_KIM_DEFAULT})
  if(DOWNLOAD_KIM)
    if(CMAKE_GENERATOR STREQUAL "Ninja")
      message(FATAL_ERROR "Cannot build downloaded KIM-API library with Ninja build tool")
    endif()
    message(STATUS "KIM-API download requested - we will build our own")
    enable_language(C)
    enable_language(Fortran)
    include(ExternalProject)
    ExternalProject_Add(kim_build
      GIT_REPOSITORY https://github.com/openkim/kim-api.git
      GIT_TAG SimulatorModels
      #URL https://s3.openkim.org/kim-api/kim-api-2.0.2.txz
      #URL_MD5 537d9c0abd30f85b875ebb584f9143fa
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
