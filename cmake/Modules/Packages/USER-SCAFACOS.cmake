if(PKG_USER-SCAFACOS)
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
    include(ExternalProject)
    ExternalProject_Add(scafacos_build
      URL https://github.com/scafacos/scafacos/releases/download/v1.0.1/scafacos-1.0.1.tar.gz
      URL_MD5 bd46d74e3296bd8a444d731bb10c1738
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --disable-doc
                                               --enable-fcs-solvers=fmm,p2nfft,direct,ewald,p3m
                                               --with-internal-fftw --with-internal-pfft
                                               --with-internal-pnfft ${CONFIGURE_REQUEST_PIC}
                                               FC=${CMAKE_MPI_Fortran_COMPILER}
                                               CXX=${CMAKE_MPI_CXX_COMPILER}
                                               CC=${CMAKE_MPI_C_COMPILER}
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
    set(SCAFACOS_BUILD_DIR ${INSTALL_DIR})
    set(SCAFACOS_INCLUDE_DIRS ${SCAFACOS_BUILD_DIR}/include)
    list(APPEND LAMMPS_DEPS scafacos_build)
    # list and order from pkg_config file of ScaFaCoS build
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_direct.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_ewald.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_fmm.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_p2nfft.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_p3m.a)
    list(APPEND LAMMPS_LINK_LIBS ${GSL_LIBRARIES})
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_near.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_gridsort.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_resort.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_redist.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_common.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_pnfft.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_pfft.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_fftw3_mpi.a)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_BUILD_DIR}/lib/libfcs_fftw3.a)
    list(APPEND LAMMPS_LINK_LIBS ${MPI_Fortran_LIBRARIES})
    list(APPEND LAMMPS_LINK_LIBS ${MPI_C_LIBRARIES})
  else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(SCAFACOS REQUIRED scafacos)
    list(APPEND LAMMPS_LINK_LIBS ${SCAFACOS_LDFLAGS})
  endif()
  include_directories(${SCAFACOS_INCLUDE_DIRS})
endif()
