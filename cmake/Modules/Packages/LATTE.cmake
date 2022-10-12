enable_language(Fortran)

# using lammps in a super-build setting
if(TARGET LATTE::latte)
  target_link_libraries(lammps PRIVATE LATTE::latte)
  return()
endif()

find_package(LATTE 1.2.2 CONFIG)
if(LATTE_FOUND)
  set(DOWNLOAD_LATTE_DEFAULT OFF)
else()
  set(DOWNLOAD_LATTE_DEFAULT ON)
endif()
option(DOWNLOAD_LATTE "Download the LATTE library instead of using an already installed one" ${DOWNLOAD_LATTE_DEFAULT})
if(DOWNLOAD_LATTE)
  message(STATUS "LATTE download requested - we will build our own")
  set(LATTE_URL "https://github.com/lanl/LATTE/archive/v1.2.2.tar.gz" CACHE STRING "URL for LATTE tarball")
  set(LATTE_MD5 "820e73a457ced178c08c71389a385de7" CACHE STRING "MD5 checksum of LATTE tarball")
  mark_as_advanced(LATTE_URL)
  mark_as_advanced(LATTE_MD5)

  # CMake cannot pass BLAS or LAPACK library variable to external project if they are a list
  list(LENGTH BLAS_LIBRARIES} NUM_BLAS)
  list(LENGTH LAPACK_LIBRARIES NUM_LAPACK)
  if((NUM_BLAS GREATER 1) OR (NUM_LAPACK GREATER 1) AND NOT USE_INTERNAL_LINALG)
    message(FATAL_ERROR "Cannot compile downloaded LATTE library due to a technical limitation. "
                        "Try to configure LAMMPS with '-D USE_INTERNAL_LINALG=on' added as a workaround.")
  endif()

  include(ExternalProject)
  ExternalProject_Add(latte_build
    URL     ${LATTE_URL}
    URL_MD5 ${LATTE_MD5}
    SOURCE_SUBDIR cmake
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ${CMAKE_REQUEST_PIC} -DCMAKE_INSTALL_LIBDIR=lib
    -DBLAS_LIBRARIES=${BLAS_LIBRARIES} -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES}
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
    -DCMAKE_Fortran_FLAGS_${BTYPE}=${CMAKE_Fortran_FLAGS_${BTYPE}} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM} -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    BUILD_BYPRODUCTS <INSTALL_DIR>/lib/liblatte.a
  )
  ExternalProject_get_property(latte_build INSTALL_DIR)
  add_library(LAMMPS::LATTE UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::LATTE PROPERTIES
    IMPORTED_LOCATION "${INSTALL_DIR}/lib/liblatte.a"
    INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
  target_link_libraries(lammps PRIVATE LAMMPS::LATTE)
  add_dependencies(LAMMPS::LATTE latte_build)
else()
  find_package(LATTE 1.2.2 REQUIRED CONFIG)
  target_link_libraries(lammps PRIVATE LATTE::latte)
endif()
