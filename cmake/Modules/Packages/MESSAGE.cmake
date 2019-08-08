if(PKG_MESSAGE)
  option(MESSAGE_ZMQ "Use ZeroMQ in MESSAGE package" OFF)
  file(GLOB_RECURSE cslib_SOURCES ${LAMMPS_LIB_SOURCE_DIR}/message/cslib/[^.]*.F
      ${LAMMPS_LIB_SOURCE_DIR}/message/cslib/[^.]*.c
      ${LAMMPS_LIB_SOURCE_DIR}/message/cslib/[^.]*.cpp)

  add_library(cslib STATIC ${cslib_SOURCES})
  if(BUILD_MPI)
    target_compile_definitions(cslib PRIVATE -DMPI_YES)
    set_target_properties(cslib PROPERTIES OUTPUT_NAME "csmpi")
  else()
    target_compile_definitions(cslib PRIVATE -DMPI_NO)
    target_include_directories(cslib PRIVATE ${LAMMPS_LIB_SOURCE_DIR}/message/cslib/src/STUBS_MPI)
    set_target_properties(cslib PROPERTIES OUTPUT_NAME "csnompi")
  endif()

  if(MESSAGE_ZMQ)
    target_compile_definitions(cslib PRIVATE -DZMQ_YES)
    find_package(ZMQ REQUIRED)
    target_include_directories(cslib PRIVATE ${ZMQ_INCLUDE_DIRS})
    target_link_libraries(cslib PUBLIC ${ZMQ_LIBRARIES})
  else()
    target_compile_definitions(cslib PRIVATE -DZMQ_NO)
    target_include_directories(cslib PRIVATE ${LAMMPS_LIB_SOURCE_DIR}/message/cslib/src/STUBS_ZMQ)
  endif()

  list(APPEND LAMMPS_LINK_LIBS cslib)
  include_directories(${LAMMPS_LIB_SOURCE_DIR}/message/cslib/src)
endif()
