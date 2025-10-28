# fix sgcmc may only be installed if also the EAM pair style from MANYBODY is installed
if(NOT PKG_MANYBODY)
  get_property(LAMMPS_FIX_HEADERS GLOBAL PROPERTY FIX)
  list(REMOVE_ITEM LAMMPS_FIX_HEADERS ${LAMMPS_SOURCE_DIR}/MC/fix_sgcmc.h)
  set_property(GLOBAL PROPERTY FIX "${LAMMPS_FIX_HEADERS}")
  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  list(REMOVE_ITEM LAMMPS_SOURCES ${LAMMPS_SOURCE_DIR}/MC/fix_sgcmc.cpp)
  set_property(TARGET lammps PROPERTY SOURCES "${LAMMPS_SOURCES}")
endif()

# fix hmc may only be installed if also fix rigid/small from RIGID is installed
if(NOT PKG_RIGID)
  get_property(LAMMPS_FIX_HEADERS GLOBAL PROPERTY FIX)
  list(REMOVE_ITEM LAMMPS_FIX_HEADERS ${LAMMPS_SOURCE_DIR}/MC/fix_hmc.h)
  set_property(GLOBAL PROPERTY FIX "${LAMMPS_FIX_HEADERS}")
  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  list(REMOVE_ITEM LAMMPS_SOURCES ${LAMMPS_SOURCE_DIR}/MC/fix_hmc.cpp)
  set_property(TARGET lammps PROPERTY SOURCES "${LAMMPS_SOURCES}")
endif()

# fix neighbor/swap may only be installed if also the VORONOI package is installed
if(NOT PKG_VORONOI)
  get_property(LAMMPS_FIX_HEADERS GLOBAL PROPERTY FIX)
  list(REMOVE_ITEM LAMMPS_FIX_HEADERS ${LAMMPS_SOURCE_DIR}/MC/fix_neighbor_swap.h)
  set_property(GLOBAL PROPERTY FIX "${LAMMPS_FIX_HEADERS}")
  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  list(REMOVE_ITEM LAMMPS_SOURCES ${LAMMPS_SOURCE_DIR}/MC/fix_neighbor_swap.cpp)
  set_property(TARGET lammps PROPERTY SOURCES "${LAMMPS_SOURCES}")
endif()
