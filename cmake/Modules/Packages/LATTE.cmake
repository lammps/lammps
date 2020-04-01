if(PKG_LATTE)
  enable_language(Fortran)
  find_package(LATTE)
  if(LATTE_FOUND)
    set(DOWNLOAD_LATTE_DEFAULT OFF)
  else()
    set(DOWNLOAD_LATTE_DEFAULT ON)
  endif()
  option(DOWNLOAD_LATTE "Download the LATTE library instead of using an already installed one" ${DOWNLOAD_LATTE_DEFAULT})
  if(DOWNLOAD_LATTE)
    message(STATUS "LATTE download requested - we will build our own")
    include(ExternalProject)
    ExternalProject_Add(latte_build
      URL https://github.com/lanl/LATTE/archive/v1.2.1.tar.gz
      URL_MD5 85ac414fdada2d04619c8f936344df14
      SOURCE_SUBDIR cmake
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ${CMAKE_REQUEST_PIC} -DCMAKE_INSTALL_LIBDIR=lib
      -DBLAS_LIBRARIES=${BLAS_LIBRARIES} -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
      -DCMAKE_Fortran_FLAGS_${BTYPE}=${CMAKE_Fortran_FLAGS_${BTYPE}} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM} -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
      BUILD_BYPRODUCTS <INSTALL_DIR>/lib/liblatte.a
    )
    ExternalProject_get_property(latte_build INSTALL_DIR)
    add_library(LAMMPS::LATTE UNKNOWN IMPORTED)
    set_target_properties(LAMMPS::LATTE PROPERTIES
      IMPORTED_LOCATION "${INSTALL_DIR}/lib/liblatte.a"
      INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
    target_link_libraries(lammps PRIVATE LAMMPS::LATTE)
    add_dependencies(LAMMPS::LATTE latte_build)
  else()
    find_package(LATTE)
    if(NOT LATTE_FOUND)
      message(FATAL_ERROR "LATTE library not found, help CMake to find it by setting LATTE_LIBRARY, or set DOWNLOAD_LATTE=ON to download it")
    endif()
    target_link_libraries(lammps PRIVATE LATTE::latte)
  endif()
endif()
