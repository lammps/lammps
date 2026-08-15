# compute xrd/fft may only be installed if also the FFT wrappers from KSPACE are installed
if(NOT PKG_KSPACE)
  get_property(LAMMPS_COMPUTE_HEADERS GLOBAL PROPERTY COMPUTE)
  list(REMOVE_ITEM LAMMPS_COMPUTE_HEADERS ${LAMMPS_SOURCE_DIR}/DIFFRACTION/compute_xrd_fft.h)
  set_property(GLOBAL PROPERTY COMPUTE "${LAMMPS_COMPUTE_HEADERS}")
  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  list(REMOVE_ITEM LAMMPS_SOURCES ${LAMMPS_SOURCE_DIR}/DIFFRACTION/compute_xrd_fft.cpp)
  set_property(TARGET lammps PROPERTY SOURCES "${LAMMPS_SOURCES}")
endif()
