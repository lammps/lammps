find_package(mdi QUIET)
if(${mdi_FOUND})
  set(DOWNLOAD_MDI_DEFAULT OFF)
else()
  set(DOWNLOAD_MDI_DEFAULT ON)
endif()
option(DOWNLOAD_MDI "Download and compile the MDI library instead of using an already installed one" ${DOWNLOAD_MDI_DEFAULT})

if(DOWNLOAD_MDI)
  message(STATUS "MDI download requested - we will build our own")
  set(MDI_URL "https://github.com/MolSSI-MDI/MDI_Library/archive/v1.2.9.tar.gz" CACHE STRING "URL for MDI tarball")
  set(MDI_MD5 "ddfa46d6ee15b4e59cfd527ec7212184" CACHE STRING "MD5 checksum for MDI tarball")
  mark_as_advanced(MDI_URL)
  mark_as_advanced(MDI_MD5)

  set(LAMMPS_LIB_MDI_BIN_DIR ${LAMMPS_LIB_BINARY_DIR}/mdi)

  include(ExternalProject)
  message(STATUS "Building mdi.")
  ExternalProject_Add(mdi_external
      URL     ${MDI_URL}
      URL_MD5 ${MDI_MD5}
      UPDATE_COMMAND ""
      CMAKE_ARGS ${CMAKE_REQUEST_PIC}
                 -DCMAKE_INSTALL_PREFIX=${LAMMPS_LIB_MDI_BIN_DIR}
                 -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                 -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                 -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                 -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                 -Dlanguage=C
      CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                       -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                       -DTargetOpenMP_FIND_COMPONENTS:STRING=C;CXX)

  # Link the lammps library against MDI
  target_include_directories(lammps PRIVATE ${LAMMPS_LIB_MDI_BIN_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/mdi)
  target_link_directories(lammps PRIVATE ${LAMMPS_LIB_MDI_BIN_DIR}/${CMAKE_INSTALL_LIBDIR}/mdi)
  target_link_libraries(lammps PRIVATE mdi)
  add_dependencies(lammps mdi_external)

  # Link the lammps executable against MDI
  target_include_directories(lmp PRIVATE ${LAMMPS_LIB_MDI_BIN_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/mdi)
  target_link_directories(lmp PRIVATE ${LAMMPS_LIB_MDI_BIN_DIR}/${CMAKE_INSTALL_LIBDIR}/mdi)
  target_link_libraries(lmp PRIVATE mdi)
  add_dependencies(lmp mdi_external)

else()

  find_package(mdi)
  if(NOT mdi_FOUND)
    message(FATAL_ERROR "MDI library not found. Help CMake to find it "
      "by setting mdi_LIBRARY and mdi_INCLUDE_DIR, or set DOWNLOAD_MDI=ON "
      "to download and compile it")
  endif()

  # Link the lammps library against MDI
  target_include_directories(lammps PRIVATE ${mdi_INCLUDE_DIR})
  target_link_libraries(lammps PRIVATE ${mdi_LIBRARY})

  # Link the lammps executable against MDI
  #target_include_directories(lmp PRIVATE ${mdi_INCLUDE_DIR})
  #target_link_libraries(lmp PRIVATE ${mdi_LIBRARY}
endif()

target_compile_definitions(lammps PRIVATE -DLMP_USER_MDI)
target_compile_definitions(lmp PRIVATE -DLMP_USER_MDI)
