if(PKG_USER-QMMM)
  enable_language(Fortran)
  enable_language(C)

  message(WARNING "Building QMMM with CMake is still experimental")
  find_package(QE REQUIRED)
  include_directories(${QE_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS ${QE_LIBRARIES})
endif()
