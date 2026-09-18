# Fix uvt uses the combined temperature compute from EXTRA-COMPUTE.
if(NOT PKG_EXTRA-COMPUTE)
  get_property(LAMMPS_FIX_HEADERS GLOBAL PROPERTY FIX)
  list(REMOVE_ITEM LAMMPS_FIX_HEADERS ${LAMMPS_SOURCE_DIR}/EXTRA-FIX/fix_uvt.h)
  set_property(GLOBAL PROPERTY FIX "${LAMMPS_FIX_HEADERS}")
  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  list(REMOVE_ITEM LAMMPS_SOURCES ${LAMMPS_SOURCE_DIR}/EXTRA-FIX/fix_uvt.cpp)
  set_property(TARGET lammps PROPERTY SOURCES "${LAMMPS_SOURCES}")
endif()
