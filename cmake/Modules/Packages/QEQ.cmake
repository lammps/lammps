# Fix qeq/fire requires MANYBODY (i.e. COMB and COMB3) to be installed
if(PKG_QEQ)
  set(QEQ_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/QEQ)
  file(GLOB QEQ_HEADERS ${QEQ_SOURCES_DIR}/fix*.h)
  file(GLOB QEQ_SOURCES ${QEQ_SOURCES_DIR}/fix*.cpp)

  if(NOT PKG_MANYBODY)
    list(REMOVE_ITEM QEQ_HEADERS ${QEQ_SOURCES_DIR}/fix_qeq_fire.h)
    list(REMOVE_ITEM QEQ_SOURCES ${QEQ_SOURCES_DIR}/fix_qeq_fire.cpp)
  endif()
  set_property(GLOBAL PROPERTY "QEQ_SOURCES" "${QEQ_SOURCES}")

  foreach(MY_HEADER ${QEQ_HEADERS})
    AddStyleHeader(${MY_HEADER} FIX)
  endforeach()

  get_property(QEQ_SOURCES GLOBAL PROPERTY QEQ_SOURCES)
  list(APPEND LIB_SOURCES ${QEQ_SOURCES})
  include_directories(${QEQ_SOURCES_DIR})
endif()
