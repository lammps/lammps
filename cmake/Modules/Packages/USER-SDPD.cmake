# Fix rigid/meso requires RIGID to be installed
if(PKG_USER-SDPD)
  set(USER-SDPD_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/USER-SDPD)

  get_property(hlist GLOBAL PROPERTY FIX)
  if(NOT PKG_RIGID)
    list(REMOVE_ITEM hlist ${USER-SDPD_SOURCES_DIR}/fix_rigid_meso.h)
    list(REMOVE_ITEM LIB_SOURCES ${USER-SDPD_SOURCES_DIR}/fix_rigid_meso.cpp)
  endif()
  set_property(GLOBAL PROPERTY FIX "${hlist}")

  include_directories(${USER-SDPD_SOURCES_DIR})
endif()
