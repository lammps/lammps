# Fix qeq/fire requires MANYBODY (i.e. COMB and COMB3) to be installed
set(QEQ_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/QEQ)

get_property(hlist GLOBAL PROPERTY FIX)
if(NOT PKG_MANYBODY)
  list(REMOVE_ITEM hlist ${QEQ_SOURCES_DIR}/fix_qeq_fire.h)
  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  list(REMOVE_ITEM LAMMPS_SOURCES ${QEQ_SOURCES_DIR}/fix_qeq_fire.cpp)
  set_property(TARGET lammps PROPERTY SOURCES ${LAMMPS_SOURCES})
endif()
set_property(GLOBAL PROPERTY FIX "${hlist}")

target_include_directories(lammps PRIVATE ${QEQ_SOURCES_DIR})
