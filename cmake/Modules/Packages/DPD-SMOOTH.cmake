# Fix rigid/meso requires RIGID to be installed
set(DPD-SMOOTH_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/DPD-SMOOTH)

get_property(hlist GLOBAL PROPERTY FIX)
if(NOT PKG_RIGID)
  list(REMOVE_ITEM hlist ${DPD-SMOOTH_SOURCES_DIR}/fix_rigid_meso.h)
  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  list(REMOVE_ITEM LAMMPS_SOURCES ${DPD-SMOOTH_SOURCES_DIR}/fix_rigid_meso.cpp)
  set_property(TARGET lammps PROPERTY SOURCES ${LAMMPS_SOURCES})
endif()
set_property(GLOBAL PROPERTY FIX "${hlist}")

target_include_directories(lammps PRIVATE ${DPD-SMOOTH_SOURCES_DIR})
