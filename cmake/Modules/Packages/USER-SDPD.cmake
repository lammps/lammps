# Fix rigid/meso requires RIGID to be installed
if(PKG_USER-SDPD)
  set(USER-SDPD_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/USER-SDPD)

  get_property(hlist GLOBAL PROPERTY FIX)
  if(NOT PKG_RIGID)
    list(REMOVE_ITEM hlist ${USER-SDPD_SOURCES_DIR}/fix_rigid_meso.h)
    get_target_property(LAMMPS_SOURCES lammps SOURCES)
    list(REMOVE_ITEM LAMMPS_SOURCES ${USER-SDPD_SOURCES_DIR}/fix_rigid_meso.cpp)
    set_property(TARGET lammps PROPERTY SOURCES ${LAMMPS_SOURCES})
  endif()
  set_property(GLOBAL PROPERTY FIX "${hlist}")

  include_directories(${USER-SDPD_SOURCES_DIR})
endif()
