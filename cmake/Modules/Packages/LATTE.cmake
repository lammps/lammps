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
    # Workaround for cross compilation with MinGW where ${CMAKE_INSTALL_LIBDIR}
    # is a full path, so we need to remove the prefix
    string(REPLACE ${CMAKE_INSTALL_PREFIX} "" _LATTE_LIBDIR ${CMAKE_INSTALL_LIBDIR})
    include(ExternalProject)
    ExternalProject_Add(latte_build
      URL https://github.com/lanl/LATTE/archive/v1.2.1.tar.gz
      URL_MD5 85ac414fdada2d04619c8f936344df14
      SOURCE_SUBDIR cmake
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ${CMAKE_REQUEST_PIC}
      -DBLAS_LIBRARIES=${BLAS_LIBRARIES} -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
      -DCMAKE_Fortran_FLAGS_${BTYPE}=${CMAKE_Fortran_FLAGS_${BTYPE}} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM} -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
      BUILD_BYPRODUCTS <INSTALL_DIR>/${_LATTE_LIBDIR}/liblatte.a
    )
    list(APPEND LAMMPS_DEPS latte_build)
    ExternalProject_get_property(latte_build INSTALL_DIR)
    set(LATTE_LIBRARIES ${INSTALL_DIR}/${_LATTE_LIBDIR}/liblatte.a)
  else()
    find_package(LATTE)
    if(NOT LATTE_FOUND)
      message(FATAL_ERROR "LATTE library not found, help CMake to find it by setting LATTE_LIBRARY, or set DOWNLOAD_LATTE=ON to download it")
    endif()
  endif()
  if(NOT LAPACK_FOUND)
    add_dependencies(latte_build linalg)
  endif()
  target_link_libraries(lammps PRIVATE ${LATTE_LIBRARIES} ${LAPACK_LIBRARIES})
endif()
