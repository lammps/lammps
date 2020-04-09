find_package(GSL REQUIRED)
find_package(MSCG QUIET)
if(MSGC_FOUND)
  set(DOWNLOAD_MSCG_DEFAULT OFF)
else()
  set(DOWNLOAD_MSCG_DEFAULT ON)
endif()
option(DOWNLOAD_MSCG "Download MSCG library instead of using an already installed one)" ${DOWNLOAD_MSCG_DEFAULT})
if(DOWNLOAD_MSCG)
  include(ExternalProject)
  ExternalProject_Add(mscg_build
    URL https://github.com/uchicago-voth/MSCG-release/archive/1.7.3.1.tar.gz
    URL_MD5 8c45e269ee13f60b303edd7823866a91
    SOURCE_SUBDIR src/CMake
    CMAKE_ARGS ${CMAKE_REQUEST_PIC} ${EXTRA_MSCG_OPTS}
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
               -DBLAS_LIBRARIES=${BLAS_LIBRARIES} -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES}
               -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
               -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    BUILD_COMMAND ${CMAKE_COMMAND} --build . --target mscg
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS <BINARY_DIR>/libmscg.a
    )
  ExternalProject_get_property(mscg_build BINARY_DIR)
  ExternalProject_get_property(mscg_build SOURCE_DIR)
  file(MAKE_DIRECTORY ${SOURCE_DIR}/src)
  add_library(LAMMPS::MSCG UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::MSCG PROPERTIES
    IMPORTED_LOCATION "${BINARY_DIR}/libmscg.a"
    INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/src"
    INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
  target_link_libraries(lammps PRIVATE LAMMPS::MSCG)
  add_dependencies(LAMMPS::MSCG mscg_build)
  if(BUILD_LIB AND NOT BUILD_SHARED_LIBS)
    install(CODE "MESSAGE(FATAL_ERROR \"Installing liblammps with downloaded libraries is currently not supported.\")")
  endif()
else()
  find_package(MSCG)
  if(NOT MSCG_FOUND)
    message(FATAL_ERROR "MSCG not found, help CMake to find it by setting MSCG_LIBRARY and MSCG_INCLUDE_DIRS, or set DOWNLOAD_MSCG=ON to download it")
  endif()
  target_link_libraries(lammps PRIVATE MSCG::MSCG)
endif()
target_link_libraries(lammps PRIVATE GSL::gsl ${LAPACK_LIBRARIES})
