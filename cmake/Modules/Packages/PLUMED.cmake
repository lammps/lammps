# Plumed2 support for PLUMED package

if(BUILD_MPI)
  set(PLUMED_CONFIG_MPI "--enable-mpi")
  set(PLUMED_CONFIG_CC  ${CMAKE_MPI_C_COMPILER})
  set(PLUMED_CONFIG_CXX  ${CMAKE_MPI_CXX_COMPILER})
  set(PLUMED_CONFIG_CPP "-I ${MPI_CXX_INCLUDE_PATH}")
  set(PLUMED_CONFIG_LIB "${MPI_CXX_LIBRARIES}")
  set(PLUMED_CONFIG_DEP "mpi4win_build")
else()
  set(PLUMED_CONFIG_MPI "--disable-mpi")
  set(PLUMED_CONFIG_CC  ${CMAKE_C_COMPILER})
  set(PLUMED_CONFIG_CXX  ${CMAKE_CXX_COMPILER})
  set(PLUMED_CONFIG_CPP "")
  set(PLUMED_CONFIG_LIB "")
  set(PLUMED_CONFIG_DEP "")
endif()
if(BUILD_OMP)
  set(PLUMED_CONFIG_OMP "--enable-openmp")
else()
  set(PLUMED_CONFIG_OMP "--disable-openmp")
endif()

set(PLUMED_URL "https://github.com/plumed/plumed2/releases/download/v2.8.3/plumed-src-2.8.3.tgz"
  CACHE STRING "URL for PLUMED tarball")
set(PLUMED_MD5 "76d23cd394eba9e6530316ed1184e219" CACHE STRING "MD5 checksum of PLUMED tarball")

mark_as_advanced(PLUMED_URL)
mark_as_advanced(PLUMED_MD5)
GetFallbackURL(PLUMED_URL PLUMED_FALLBACK)

if((CMAKE_SYSTEM_NAME STREQUAL "Windows") AND (CMAKE_CROSSCOMPILING))
  if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    set(CROSS_CONFIGURE mingw64-configure)
  elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86")
    set(CROSS_CONFIGURE mingw32-configure)
  else()
    message(FATAL_ERROR "Unsupported target system: ${CMAKE_SYSTEM_NAME}/${CMAKE_SYSTEM_PROCESSOR}")
  endif()
  message(STATUS "Downloading and cross-compiling Plumed2 for ${CMAKE_SYSTEM_NAME}/${CMAKE_SYSTEM_PROCESSOR} with ${CROSS_CONFIGURE}")
  include(ExternalProject)
  ExternalProject_Add(plumed_build
    URL     ${PLUMED_URL} ${PLUMED_FALLBACK}
    URL_MD5 ${PLUMED_MD5}
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ${CROSS_CONFIGURE} --disable-shared --disable-bsymbolic
                                         --disable-python --enable-cxx=11
                                         --enable-modules=-adjmat:+crystallization:-dimred:+drr:+eds:-fisst:+funnel:+logmfd:+manyrestraints:+maze:+opes:+multicolvar:-pamm:-piv:+s2cm:-sasa:-ves
                                         ${PLUMED_CONFIG_OMP}
                                         ${PLUMED_CONFIG_MPI}
                                         CXX=${PLUMED_CONFIG_CXX}
                                         CC=${PLUMED_CONFIG_CC}
                                         CPPFLAGS=${PLUMED_CONFIG_CPP}
                                         LIBS=${PLUMED_CONFIG_LIB}
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS "<SOURCE_DIR>/src/lib/install/libplumed.a" "<SOURCE_DIR>/src/lib/install/plumed.exe"
    DEPENDS "${PLUMED_MPI_CONFIG_DEP}"
  )
  ExternalProject_Get_Property(plumed_build SOURCE_DIR)
  set(PLUMED_BUILD_DIR ${SOURCE_DIR})
  set(PLUMED_INSTALL_DIR ${PLUMED_BUILD_DIR}/src/lib/install)
  file(MAKE_DIRECTORY ${PLUMED_BUILD_DIR}/src/include)

  add_library(LAMMPS::PLUMED UNKNOWN IMPORTED)
  add_dependencies(LAMMPS::PLUMED plumed_build)
  set_target_properties(LAMMPS::PLUMED PROPERTIES
    IMPORTED_LOCATION "${PLUMED_INSTALL_DIR}/libplumed.a"
    INTERFACE_LINK_LIBRARIES "-Wl,--image-base -Wl,0x10000000 -lfftw3 -lz  -fstack-protector -lssp -fopenmp"
    INTERFACE_INCLUDE_DIRECTORIES "${PLUMED_BUILD_DIR}/src/include")

  add_custom_target(plumed_copy ALL ${CMAKE_COMMAND} -E rm -rf ${CMAKE_BINARY_DIR}/plumed.exe ${CMAKE_BINARY_DIR}/plumed_patches
    COMMAND ${CMAKE_COMMAND} -E copy ${PLUMED_INSTALL_DIR}/plumed.exe ${CMAKE_BINARY_DIR}/plumed.exe
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PLUMED_BUILD_DIR}/patches ${CMAKE_BINARY_DIR}/patches
    BYPRODUCTS ${CMAKE_BINARY_DIR}/plumed.exe ${CMAKE_BINARY_DIR}/patches
    DEPENDS plumed_build
    COMMENT "Copying Plumed files"
  )

else()

  set(PLUMED_MODE "static" CACHE STRING "Linkage mode for Plumed2 library")
  set(PLUMED_MODE_VALUES static shared runtime)
  set_property(CACHE PLUMED_MODE PROPERTY STRINGS ${PLUMED_MODE_VALUES})
  validate_option(PLUMED_MODE PLUMED_MODE_VALUES)
  string(TOUPPER ${PLUMED_MODE} PLUMED_MODE)

  set(PLUMED_LINK_LIBS)
  if(PLUMED_MODE STREQUAL "STATIC")
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
    find_package(GSL REQUIRED)
    list(APPEND PLUMED_LINK_LIBS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} GSL::gsl)
    find_package(ZLIB QUIET)
    if(ZLIB_FOUND)
      list(APPEND PLUMED_LINK_LIBS ZLIB::ZLIB)
    endif()
    find_package(FFTW3 QUIET)
    if(FFTW3_FOUND)
      list(APPEND PLUMED_LINK_LIBS FFTW3::FFTW3)
    endif()
  endif()

  find_package(PkgConfig QUIET)
  set(DOWNLOAD_PLUMED_DEFAULT ON)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(PLUMED QUIET plumed)
    if(PLUMED_FOUND)
      set(DOWNLOAD_PLUMED_DEFAULT OFF)
    endif()
  endif()

  option(DOWNLOAD_PLUMED "Download Plumed package instead of using an already installed one" ${DOWNLOAD_PLUMED_DEFAULT})
  if(DOWNLOAD_PLUMED)
    message(STATUS "PLUMED download requested - we will build our own")
    if(PLUMED_MODE STREQUAL "STATIC")
      set(PLUMED_BUILD_BYPRODUCTS "<INSTALL_DIR>/lib/${CMAKE_STATIC_LIBRARY_PREFIX}plumed${CMAKE_STATIC_LIBRARY_SUFFIX}")
    elseif(PLUMED_MODE STREQUAL "SHARED")
      set(PLUMED_BUILD_BYPRODUCTS "<INSTALL_DIR>/lib/${CMAKE_SHARED_LIBRARY_PREFIX}plumed${CMAKE_SHARED_LIBRARY_SUFFIX};<INSTALL_DIR>/lib/${CMAKE_SHARED_LIBRARY_PREFIX}plumedKernel${CMAKE_SHARED_LIBRARY_SUFFIX}")
    elseif(PLUMED_MODE STREQUAL "RUNTIME")
      set(PLUMED_BUILD_BYPRODUCTS "<INSTALL_DIR>/lib/${CMAKE_STATIC_LIBRARY_PREFIX}plumedWrapper${CMAKE_STATIC_LIBRARY_PREFIX}")
    endif()

    include(ExternalProject)
    ExternalProject_Add(plumed_build
      URL     ${PLUMED_URL} ${PLUMED_FALLBACK}
      URL_MD5 ${PLUMED_MD5}
      BUILD_IN_SOURCE 1
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
                                             ${CONFIGURE_REQUEST_PIC}
                                             --enable-modules=all
                                             --enable-cxx=11
                                             --disable-python
                                             ${PLUMED_CONFIG_MPI}
                                             ${PLUMED_CONFIG_OMP}
                                             CXX=${PLUMED_CONFIG_CXX}
                                             CC=${PLUMED_CONFIG_CC}
      BUILD_BYPRODUCTS ${PLUMED_BUILD_BYPRODUCTS}
    )
    ExternalProject_get_property(plumed_build INSTALL_DIR)
    add_library(LAMMPS::PLUMED UNKNOWN IMPORTED)
    add_dependencies(LAMMPS::PLUMED plumed_build)
    if(PLUMED_MODE STREQUAL "STATIC")
      set_target_properties(LAMMPS::PLUMED PROPERTIES IMPORTED_LOCATION ${INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}plumed${CMAKE_STATIC_LIBRARY_SUFFIX} INTERFACE_LINK_LIBRARIES "${PLUMED_LINK_LIBS};${CMAKE_DL_LIBS}")
    elseif(PLUMED_MODE STREQUAL "SHARED")
      set_target_properties(LAMMPS::PLUMED PROPERTIES IMPORTED_LOCATION ${INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}plumed${CMAKE_SHARED_LIBRARY_SUFFIX} INTERFACE_LINK_LIBRARIES "${INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}plumedKernel${CMAKE_SHARED_LIBRARY_SUFFIX};${CMAKE_DL_LIBS}")
    elseif(PLUMED_MODE STREQUAL "RUNTIME")
      set_target_properties(LAMMPS::PLUMED PROPERTIES INTERFACE_COMPILE_DEFINITIONS "__PLUMED_DEFAULT_KERNEL=${INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}plumedKernel${CMAKE_SHARED_LIBRARY_SUFFIX}")
      set_target_properties(LAMMPS::PLUMED PROPERTIES IMPORTED_LOCATION ${INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}plumedWrapper${CMAKE_STATIC_LIBRARY_SUFFIX} INTERFACE_LINK_LIBRARIES "${CMAKE_DL_LIBS}")
    endif()
    set_target_properties(LAMMPS::PLUMED PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_DIR}/include)
    file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
  else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(PLUMED REQUIRED plumed)
    add_library(LAMMPS::PLUMED INTERFACE IMPORTED)
    if(PLUMED_MODE STREQUAL "STATIC")
      include(${PLUMED_LIBDIR}/plumed/src/lib/Plumed.cmake.static)
    elseif(PLUMED_MODE STREQUAL "SHARED")
      include(${PLUMED_LIBDIR}/plumed/src/lib/Plumed.cmake.shared)
    elseif(PLUMED_MODE STREQUAL "RUNTIME")
      set_target_properties(LAMMPS::PLUMED PROPERTIES INTERFACE_COMPILE_DEFINITIONS "__PLUMED_DEFAULT_KERNEL=${PLUMED_LIBDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}plumedKernel${CMAKE_SHARED_LIBRARY_SUFFIX}")
      include(${PLUMED_LIBDIR}/plumed/src/lib/Plumed.cmake.runtime)
    endif()
    set_target_properties(LAMMPS::PLUMED PROPERTIES INTERFACE_LINK_LIBRARIES "${PLUMED_LOAD}")
    set_target_properties(LAMMPS::PLUMED PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${PLUMED_INCLUDE_DIRS}")
  endif()
endif()

target_link_libraries(lammps PRIVATE LAMMPS::PLUMED)
