if(PKG_USER-MOLFILE)
  if (CMAKE_VERSION VERSION_LESS "3.10") # due to INTERFACE without a library
    message(FATAL_ERROR "For configuring USER-MOLFILE you need CMake 3.10 or later")
  endif()

  set(MOLFILE_INCLUDE_DIRS "${LAMMPS_LIB_SOURCE_DIR}/molfile" CACHE STRING "Path to VMD molfile plugin headers")
  add_library(molfile INTERFACE)
  target_include_directories(molfile INTERFACE ${MOLFILE_INCLUDE_DIRS})
  # no need to link with -ldl on windows
  if(NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    target_link_libraries(molfile INTERFACE ${CMAKE_DL_LIBS})
  endif()
  list(APPEND LAMMPS_LINK_LIBS molfile)
endif()
