# pair style dpd/coul/slater/long may only be installed if also KSPACE is installed
if(NOT PKG_KSPACE)
  get_property(LAMMPS_PAIR_HEADERS GLOBAL PROPERTY PAIR)
  list(REMOVE_ITEM LAMMPS_PAIR_HEADERS ${LAMMPS_SOURCE_DIR}/DPD-BASIC/pair_dpd_coul_slater_long.h)
  set_property(GLOBAL PROPERTY PAIR "${LAMMPS_PAIR_HEADERS}")
  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  list(REMOVE_ITEM LAMMPS_SOURCES ${LAMMPS_SOURCE_DIR}/DPD-BASIC/pair_dpd_coul_slater_long.cpp)
  set_property(TARGET lammps PROPERTY SOURCES "${LAMMPS_SOURCES}")
endif()
