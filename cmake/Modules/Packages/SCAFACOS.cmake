enable_language(Fortran)
enable_language(C)

find_package(GSL REQUIRED)
find_package(PkgConfig QUIET)
find_package(MPI REQUIRED COMPONENTS C MPICXX CXX Fortran)
set(DOWNLOAD_SCAFACOS_DEFAULT ON)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(SCAFACOS QUIET scafacos)
  if(SCAFACOS_FOUND)
    set(DOWNLOAD_SCAFACOS_DEFAULT OFF)
  endif()
endif()
option(DOWNLOAD_SCAFACOS "Download ScaFaCoS library instead of using an already installed one" ${DOWNLOAD_SCAFACOS_DEFAULT})
if(DOWNLOAD_SCAFACOS)
  message(STATUS "ScaFaCoS download requested - we will build our own")
  set(SCAFACOS_URL "https://github.com/scafacos/scafacos/releases/download/v1.0.4/scafacos-1.0.4.tar.gz" CACHE STRING "URL for SCAFACOS tarball")
  set(SCAFACOS_SHA256 "6634c4202e825e771d1dd75bbe9cac5cee41136c87653fde98fbd634681c1be6" CACHE STRING "SHA256 checksum of SCAFACOS tarball")
  mark_as_advanced(SCAFACOS_URL)
  mark_as_advanced(SCAFACOS_SHA256)
  GetFallbackURL(SCAFACOS_URL SCAFACOS_FALLBACK)
  set(SCAFACOS_CXX_FLAGS "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}} ${CMAKE_CXX_FLAGS}")
  set(SCAFACOS_C_FLAGS "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}} ${CMAKE_C_FLAGS}")

  include(ExternalProject)
  ExternalProject_Add(scafacos_build
    URL     ${SCAFACOS_URL} ${SCAFACOS_FALLBACK}
    URL_HASH SHA256=${SCAFACOS_SHA256}
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --disable-doc
                                             --enable-fcs-solvers=fmm,p2nfft,direct,ewald,p3m
                                             --with-internal-fftw --with-internal-pfft
                                             --with-internal-pnfft ${CONFIGURE_REQUEST_PIC}
                                             FC=${CMAKE_MPI_Fortran_COMPILER}
                                             CXX=${CMAKE_MPI_CXX_COMPILER}
                                             CC=${CMAKE_MPI_C_COMPILER}
                                             F77=
                                             CFLAGS=${SCAFACOS_C_FLAGS}
                                             CXXFLAGS=${SCAFACOS_CXX_FLAGS}
    BUILD_BYPRODUCTS
      <INSTALL_DIR>/lib/libfcs.a
      <INSTALL_DIR>/lib/libfcs_direct.a
      <INSTALL_DIR>/lib/libfcs_ewald.a
      <INSTALL_DIR>/lib/libfcs_fmm.a
      <INSTALL_DIR>/lib/libfcs_p2nfft.a
      <INSTALL_DIR>/lib/libfcs_p3m.a
      <INSTALL_DIR>/lib/libfcs_near.a
      <INSTALL_DIR>/lib/libfcs_gridsort.a
      <INSTALL_DIR>/lib/libfcs_resort.a
      <INSTALL_DIR>/lib/libfcs_redist.a
      <INSTALL_DIR>/lib/libfcs_common.a
      <INSTALL_DIR>/lib/libfcs_pnfft.a
      <INSTALL_DIR>/lib/libfcs_pfft.a
      <INSTALL_DIR>/lib/libfcs_fftw3_mpi.a
      <INSTALL_DIR>/lib/libfcs_fftw3.a
  )
  ExternalProject_get_property(scafacos_build INSTALL_DIR)
  file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
  add_library(LAMMPS::SCAFACOS UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::SCAFACOS PROPERTIES
    IMPORTED_LOCATION "${INSTALL_DIR}/lib/libfcs.a"
    INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/include"
    INTERFACE_LINK_LIBRARIES "${INSTALL_DIR}/lib/libfcs.a;${INSTALL_DIR}/lib/libfcs_direct.a;${INSTALL_DIR}/lib/libfcs_ewald.a;${INSTALL_DIR}/lib/libfcs_fmm.a;${INSTALL_DIR}/lib/libfcs_p2nfft.a;${INSTALL_DIR}/lib/libfcs_p3m.a;GSL::gsl;${INSTALL_DIR}/lib/libfcs_near.a;${INSTALL_DIR}/lib/libfcs_gridsort.a;${INSTALL_DIR}/lib/libfcs_resort.a;${INSTALL_DIR}/lib/libfcs_redist.a;${INSTALL_DIR}/lib/libfcs_common.a;${INSTALL_DIR}/lib/libfcs_pnfft.a;${INSTALL_DIR}/lib/libfcs_pfft.a;${INSTALL_DIR}/lib/libfcs_fftw3_mpi.a;${INSTALL_DIR}/lib/libfcs_fftw3.a;MPI::MPI_CXX;MPI::MPI_Fortran;MPI::MPI_C;${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
  target_link_libraries(lammps PRIVATE LAMMPS::SCAFACOS)
  add_dependencies(LAMMPS::SCAFACOS scafacos_build)
else()
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(SCAFACOS REQUIRED IMPORTED_TARGET scafacos)
  target_link_libraries(lammps PRIVATE PkgConfig::SCAFACOS)
endif()
