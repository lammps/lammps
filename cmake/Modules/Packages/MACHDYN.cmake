
# set policy to use the time of extraction as timestamps of files unpacked from downloaded
# archives, so that updating an archive version triggers rebuilding all dependent objects
 if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

find_package(Eigen3 NO_MODULE)
if(EIGEN3_FOUND)
  set(DOWNLOAD_EIGEN3_DEFAULT OFF)
else()
  set(DOWNLOAD_EIGEN3_DEFAULT ON)
endif()
option(DOWNLOAD_EIGEN3 "Download Eigen3 instead of using an already installed one)" ${DOWNLOAD_EIGEN3_DEFAULT})
if(DOWNLOAD_EIGEN3)
  message(STATUS "Eigen3 download requested - we will build our own")

  SetDownloadSettings(EIGEN3 "Eigen3"
    "${LAMMPS_THIRDPARTY_URL}/eigen-3.4.0.tar.gz"
    "8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72")
  include(ExternalProject)
  ExternalProject_Add(Eigen3_build
    URL     ${EIGEN3_URL}
    URL_HASH SHA256=${EIGEN3_SHA256}
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND ""
  )
  ExternalProject_get_property(Eigen3_build SOURCE_DIR)
  add_library(LAMMPS::EIGEN3 INTERFACE IMPORTED)
  set_target_properties(LAMMPS::EIGEN3 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}")
  target_link_libraries(lammps PRIVATE LAMMPS::EIGEN3)
  add_dependencies(LAMMPS::EIGEN3 Eigen3_build)
else()
  find_package(Eigen3 NO_MODULE)
  mark_as_advanced(Eigen3_DIR)
  if(NOT EIGEN3_FOUND)
    message(FATAL_ERROR "Eigen3 not found, help CMake to find it by setting EIGEN3_INCLUDE_DIR, or set DOWNLOAD_EIGEN3=ON to download it")
  endif()
  target_link_libraries(lammps PRIVATE Eigen3::Eigen)
endif()

# PGI/Nvidia compiler internals collide with vector intrinsics support in Eigen3
if((CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR (CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC"))
  target_compile_definitions(lammps PRIVATE -DEIGEN_DONT_VECTORIZE)
endif()

target_compile_definitions(lammps PRIVATE -DEIGEN_NO_CUDA)
