find_package(mdi)
if(${mdi_FOUND})
  set(DOWNLOAD_MDI_DEFAULT OFF)
else()
  set(DOWNLOAD_MDI_DEFAULT ON)
endif()
option(DOWNLOAD_MDI "Download and compile the MDI library instead of using an already installed one" ${DOWNLOAD_MDI_DEFAULT})
if(DOWNLOAD_MDI)
  message(STATUS "MDI download requested - we will build our own")
  set(mdi_URL "https://github.com/MolSSI-MDI/MDI_Library/archive/v1.2.7.tar.gz" CACHE STRING "URL for MDI tarball")
  set(mdi_MD5 "2f3177b30ccdbd6ae28ea3bdd5fed0db" CACHE STRING "MD5 checksum for MDI tarball")
  mark_as_advanced(mdi_URL)
  mark_as_advanced(mdi_MD5)

  set(LAMMPS_LIB_MDI_BIN_DIR ${LAMMPS_LIB_BINARY_DIR}/mdi)

  include(ExternalProject)
  message(STATUS "Building mdi.")
  ExternalProject_Add(mdi_external
      URL     ${mdi_URL}
      URL_MD5 ${mdi_MD5}
      UPDATE_COMMAND ""
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${LAMMPS_LIB_MDI_BIN_DIR}
                 -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                 -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                 -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                 -DENABLE_OPENMP=${ENABLE_OPENMP}
                 -DENABLE_XHOST=${ENABLE_XHOST}
                 -DBUILD_FPIC=${BUILD_FPIC}
                 -DENABLE_GENERIC=${ENABLE_GENERIC}
                 -DLIBC_INTERJECT=${LIBC_INTERJECT}
                 -Dlanguage=C
      CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                       -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                       -DTargetOpenMP_FIND_COMPONENTS:STRING=C;CXX)

  # Link the lammps library against MDI
  target_include_directories(lammps PRIVATE ${LAMMPS_LIB_MDI_BIN_DIR}/include/mdi)
  target_link_directories(lammps PRIVATE ${LAMMPS_LIB_MDI_BIN_DIR}/lib/mdi)
  target_link_libraries(lammps PRIVATE mdi)

  # Link the lammps executable against MDI
  target_include_directories(lmp PRIVATE ${LAMMPS_LIB_MDI_BIN_DIR}/include/mdi)
  target_link_directories(lmp PRIVATE ${LAMMPS_LIB_MDI_BIN_DIR}/lib/mdi)
  target_link_libraries(lmp PRIVATE mdi)

  add_definitions(-DLMP_USER_MDI=1)

else()

  find_package(mdi)
  if(NOT mdi_FOUND)
    message(FATAL_ERROR "MDI library not found. Help CMake to find it by setting mdi_LIBRARY and mdi_INCLUDE_DIR, or set DOWNLOAD_MDI=ON to download it")
  endif()

  # Link the lammps library against MDI
  target_include_directories(lammps PRIVATE ${mdi_INCLUDE_DIR})
  target_link_libraries(lammps PRIVATE ${mdi_LIBRARIES})

  # Link the lammps executable against MDI
  target_include_directories(lmp PRIVATE ${mdi_INCLUDE_DIR})
  target_link_libraries(lmp PRIVATE ${mdi_LIBRARIES})

endif()
