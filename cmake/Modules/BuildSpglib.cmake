# Download and compile spglib as a static library for the phana and find-symmetry tools.
# Defines the imported target SPGLIB::SYMSPG.  Including this file again after the target
# exists does nothing, so several tools can share one spglib build.

if(NOT TARGET SPGLIB::SYMSPG)
  set(SPGLIB_URL "https://github.com/spglib/spglib/archive/refs/tags/v1.11.2.1.tar.gz" CACHE STRING "URL for spglib v1.x tarball")
  set(SPGLIB_SHA256 "f6795523a04871e012e7f5f5ab97b249fa36657b73cdc9b4ea53ef023cfcaac4" CACHE STRING "SHA256 checksum for spglib tarball")
  mark_as_advanced(SPGLIB_URL)
  mark_as_advanced(SPGLIB_SHA256)
  include(LAMMPSUtils)
  GetFallbackURL(SPGLIB_URL SPGLIB_FALLBACK)

  set(SPGLIB_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/spglib_build_ext)
  string(REPLACE ";" "$<SEMICOLON>" CMAKE_OSX_ARCHITECTURES_ "${CMAKE_OSX_ARCHITECTURES}")
  include(ExternalProject)
  ExternalProject_Add(spglib_build
    URL ${SPGLIB_URL} ${SPGLIB_FALLBACK}
    URL_HASH SHA256=${SPGLIB_SHA256}
    PREFIX ${SPGLIB_BUILD_DIR}
    EXCLUDE_FROM_ALL TRUE
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SPGLIB_BUILD_DIR}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES_}
    -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
    UPDATE_COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_LIST_DIR}/../../tools/phonon/CMakeLists.spglib ${SPGLIB_BUILD_DIR}/src/spglib_build/CMakeLists.txt
    INSTALL_COMMAND ${CMAKE_COMMAND} --install ${SPGLIB_BUILD_DIR}/src/spglib_build-build --prefix ${SPGLIB_BUILD_DIR}
    BUILD_BYPRODUCTS "${SPGLIB_BUILD_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}symspg${CMAKE_STATIC_LIBRARY_SUFFIX}"
  )

  # workaround for older CMake versions
  file(MAKE_DIRECTORY ${SPGLIB_BUILD_DIR}/lib)
  file(MAKE_DIRECTORY ${SPGLIB_BUILD_DIR}/include)

  add_library(SPGLIB::SYMSPG UNKNOWN IMPORTED GLOBAL)
  add_dependencies(SPGLIB::SYMSPG spglib_build)
  set_target_properties(SPGLIB::SYMSPG PROPERTIES
    IMPORTED_LOCATION "${SPGLIB_BUILD_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}symspg${CMAKE_STATIC_LIBRARY_SUFFIX}"
    INTERFACE_INCLUDE_DIRECTORIES ${SPGLIB_BUILD_DIR}/include
  )
endif()
