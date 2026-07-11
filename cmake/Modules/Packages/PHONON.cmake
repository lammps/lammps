# fix phonon may only be installed if also the FFT wrappers from KSPACE are installed
if(NOT PKG_KSPACE)
  get_property(LAMMPS_FIX_HEADERS GLOBAL PROPERTY FIX)
  list(REMOVE_ITEM LAMMPS_FIX_HEADERS ${LAMMPS_SOURCE_DIR}/PHONON/fix_phonon.h)
  set_property(GLOBAL PROPERTY FIX "${LAMMPS_FIX_HEADERS}")
  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  list(REMOVE_ITEM LAMMPS_SOURCES ${LAMMPS_SOURCE_DIR}/PHONON/fix_phonon.cpp)
  set_property(TARGET lammps PROPERTY SOURCES "${LAMMPS_SOURCES}")
endif()
