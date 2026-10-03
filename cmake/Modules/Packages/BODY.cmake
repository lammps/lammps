# pair style body/rounded/polyhedron/granular may only be installed if also the GRANULAR package is installed
if(NOT PKG_GRANULAR)
  get_property(LAMMPS_PAIR_HEADERS GLOBAL PROPERTY PAIR)
  list(REMOVE_ITEM LAMMPS_PAIR_HEADERS ${LAMMPS_SOURCE_DIR}/BODY/pair_body_rounded_polyhedron_granular.h)
  set_property(GLOBAL PROPERTY PAIR "${LAMMPS_PAIR_HEADERS}")
  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  list(REMOVE_ITEM LAMMPS_SOURCES ${LAMMPS_SOURCE_DIR}/BODY/pair_body_rounded_polyhedron_granular.cpp)
  set_property(TARGET lammps PROPERTY SOURCES "${LAMMPS_SOURCES}")
endif()
