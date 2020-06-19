enable_language(Fortran)
enable_language(C)

find_package(GSL REQUIRED)
find_package(PkgConfig QUIET)
find_package(MPI REQUIRED)
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
  # create variables to pass our compiler flags along to the subsystem compile
  # need to apply -fallow-argument-mismatch, if the fortran compiler supports it
  include(CheckFortranCompilerFlag)
  check_fortran_compiler_flag("-fallow-argument-mismatch" GNUFortran_ARGUMENT_MISMATCH_FLAG)
  if(GNUFortran_ARGUMENT_MISMATCH_FLAG)
    set(APPEND_Fortran_FLAG "-fallow-argument-mismatch")
  endif()
  if(CMAKE_Fortran_FLAGS)
    set(SCAFACOS_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${APPEND_Fortran_FLAG}")
  else()
    set(SCAFACOS_Fortran_FLAGS "${CMAKE_Fortran_${CMAKE_BUILD_TYPE}_FLAGS} ${APPEND_Fortran_FLAG}")
  endif()
  if(CMAKE_CXX_FLAGS)
      set(SCAFACOS_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
      set(SCAFACOS_CXX_FLAGS "${CMAKE_CXX_${CMAKE_BUILD_TYPE}_FLAGS}")
  endif()
  if(CMAKE_C_FLAGS)
      set(SCAFACOS_C_FLAGS "${CMAKE_C_FLAGS}")
    else()
      set(SCAFACOS_C_FLAGS "${CMAKE_C_${CMAKE_BUILD_TYPE}_FLAGS}")
  endif()

  include(ExternalProject)
  ExternalProject_Add(scafacos_build
    URL https://github.com/scafacos/scafacos/releases/download/v1.0.1/scafacos-1.0.1.tar.gz
    URL_MD5 bd46d74e3296bd8a444d731bb10c1738
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --disable-doc
                                             --enable-fcs-solvers=fmm,p2nfft,direct,ewald,p3m
                                             --with-internal-fftw --with-internal-pfft
                                             --with-internal-pnfft ${CONFIGURE_REQUEST_PIC}
                                             FC=${CMAKE_MPI_Fortran_COMPILER} FCFLAGS=${SCAFACOS_Fortran_FLAGS}
                                             CXX=${CMAKE_MPI_CXX_COMPILER} CXXFLAGS=${SCAFACOS_CXX_FLAGS}
                                             CC=${CMAKE_MPI_C_COMPILER} CFLAGS=${SCAFACOS_C_FLAGS}
                                             F77=
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
    INTERFACE_LINK_LIBRARIES "${INSTALL_DIR}/lib/libfcs.a;${INSTALL_DIR}/lib/libfcs_direct.a;${INSTALL_DIR}/lib/libfcs_ewald.a;${INSTALL_DIR}/lib/libfcs_fmm.a;${INSTALL_DIR}/lib/libfcs_p2nfft.a;${INSTALL_DIR}/lib/libfcs_p3m.a;GSL::gsl;${INSTALL_DIR}/lib/libfcs_near.a;${INSTALL_DIR}/lib/libfcs_gridsort.a;${INSTALL_DIR}/lib/libfcs_resort.a;${INSTALL_DIR}/lib/libfcs_redist.a;${INSTALL_DIR}/lib/libfcs_common.a;${INSTALL_DIR}/lib/libfcs_pnfft.a;${INSTALL_DIR}/lib/libfcs_pfft.a;${INSTALL_DIR}/lib/libfcs_fftw3_mpi.a;${INSTALL_DIR}/lib/libfcs_fftw3.a;MPI::MPI_Fortran;MPI::MPI_C")
  target_link_libraries(lammps PRIVATE LAMMPS::SCAFACOS)
  add_dependencies(LAMMPS::SCAFACOS scafacos_build)
else()
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(SCAFACOS REQUIRED IMPORTED_TARGET scafacos)
  target_link_libraries(lammps PRIVATE PkgConfig::SCAFACOS)
endif()
