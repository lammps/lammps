find_package(Eigen3 NO_MODULE)
if(EIGEN3_FOUND)
  set(DOWNLOAD_EIGEN3_DEFAULT OFF)
else()
  set(DOWNLOAD_EIGEN3_DEFAULT ON)
endif()
option(DOWNLOAD_EIGEN3 "Download Eigen3 instead of using an already installed one)" ${DOWNLOAD_EIGEN3_DEFAULT})
if(DOWNLOAD_EIGEN3)
  message(STATUS "Eigen3 download requested - we will build our own")

  set(EIGEN3_URL "${LAMMPS_THIRDPARTY_URL}/eigen-3.4.0.tar.gz" CACHE STRING "URL for Eigen3 tarball")
  set(EIGEN3_MD5 "4c527a9171d71a72a9d4186e65bea559" CACHE STRING "MD5 checksum of Eigen3 tarball")
  mark_as_advanced(EIGEN3_URL)
  mark_as_advanced(EIGEN3_MD5)
  include(ExternalProject)
  ExternalProject_Add(Eigen3_build
    URL     ${EIGEN3_URL}
    URL_MD5 ${EIGEN3_MD5}
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
