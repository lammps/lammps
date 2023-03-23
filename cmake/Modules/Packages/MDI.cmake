find_package(mdi QUIET)
if(${mdi_FOUND})
  set(DOWNLOAD_MDI_DEFAULT OFF)
else()
  set(DOWNLOAD_MDI_DEFAULT ON)
endif()
option(DOWNLOAD_MDI "Download and compile the MDI library instead of using an already installed one" ${DOWNLOAD_MDI_DEFAULT})

if(DOWNLOAD_MDI)
  message(STATUS "MDI download requested - we will build our own")
  set(MDI_URL "https://github.com/MolSSI-MDI/MDI_Library/archive/v1.4.16.tar.gz" CACHE STRING "URL for MDI tarball")
  set(MDI_MD5 "407db44e2d79447ab5c1233af1965f65" CACHE STRING "MD5 checksum for MDI tarball")
  mark_as_advanced(MDI_URL)
  mark_as_advanced(MDI_MD5)
  GetFallbackURL(MDI_URL MDI_FALLBACK)
  enable_language(C)

  # only ON/OFF are allowed for "mpi" flag when building MDI library
  # so translate boolean value of BUILD_MPI
  # always disable MPI when cross-compiling to Windows.
  if((BUILD_MPI) AND NOT((CMAKE_SYSTEM_NAME STREQUAL "Windows") AND CMAKE_CROSSCOMPILING))
    set(MDI_USE_MPI ON)
  else()
    set(MDI_USE_MPI OFF)
  endif()

  # detect if we have python development support and thus can enable python plugins
  set(MDI_USE_PYTHON_PLUGINS OFF)
  if(CMAKE_VERSION VERSION_LESS 3.12)
    if(NOT PYTHON_VERSION_STRING)
      set(Python_ADDITIONAL_VERSIONS 3.12 3.11 3.10 3.9 3.8 3.7 3.6)
      # search for interpreter first, so we have a consistent library
      find_package(PythonInterp) # Deprecated since version 3.12
      if(PYTHONINTERP_FOUND)
        set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
      endif()
    endif()
    # search for the library matching the selected interpreter
    set(Python_ADDITIONAL_VERSIONS ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
    find_package(PythonLibs QUIET) # Deprecated since version 3.12
    if(PYTHONLIBS_FOUND)
      if(NOT (PYTHON_VERSION_STRING STREQUAL PYTHONLIBS_VERSION_STRING))
        message(FATAL_ERROR "Python Library version ${PYTHONLIBS_VERSION_STRING} does not match Interpreter version ${PYTHON_VERSION_STRING}")
      endif()
      set(MDI_USE_PYTHON_PLUGINS ON)
    endif()
  else()
    find_package(Python QUIET COMPONENTS Development)
    if(Python_Development_FOUND)
      set(MDI_USE_PYTHON_PLUGINS ON)
    endif()
  endif()
  # python plugins are not supported and thus must be always off on Windows
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    unset(Python_Development_FOUND)
    set(MDI_USE_PYTHON_PLUGINS OFF)
    if(CMAKE_CROSSCOMPILING)
      set(CMAKE_INSTALL_LIBDIR lib)
    endif()
  endif()

  # download/ build MDI library
  # always build static library with -fpic
  # support cross-compilation and ninja-build
  include(ExternalProject)
  ExternalProject_Add(mdi_build
    URL     ${MDI_URL} ${MDI_FALLBACK}
    URL_MD5 ${MDI_MD5}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/mdi_build_ext
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/mdi_build_ext
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    -Dlanguage=C
    -Dlibtype=STATIC
    -Dmpi=${MDI_USE_MPI}
    -Dplugins=ON
    -Dpython_plugins=${MDI_USE_PYTHON_PLUGINS}
    UPDATE_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/mdi_build_ext/src/mdi_build-build --target install
    BUILD_BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/mdi_build_ext/${CMAKE_INSTALL_LIBDIR}/mdi/${CMAKE_STATIC_LIBRARY_PREFIX}mdi${CMAKE_STATIC_LIBRARY_SUFFIX}"
    )

  # where is the compiled library?
  ExternalProject_get_property(mdi_build PREFIX)
  # workaround for older CMake versions
  file(MAKE_DIRECTORY ${PREFIX}/${CMAKE_INSTALL_LIBDIR}/mdi)
  file(MAKE_DIRECTORY ${PREFIX}/include/mdi)

  # create imported target for the MDI library
  add_library(LAMMPS::MDI UNKNOWN IMPORTED)
  add_dependencies(LAMMPS::MDI mdi_build)
  set_target_properties(LAMMPS::MDI PROPERTIES
    IMPORTED_LOCATION "${PREFIX}/${CMAKE_INSTALL_LIBDIR}/mdi/${CMAKE_STATIC_LIBRARY_PREFIX}mdi${CMAKE_STATIC_LIBRARY_SUFFIX}"
    INTERFACE_INCLUDE_DIRECTORIES ${PREFIX}/include/mdi
  )

  set(MDI_DEP_LIBS "")
  # if compiling with python plugins we need
  # to add python libraries as dependency.
  if(MDI_USE_PYTHON_PLUGINS)
    if(CMAKE_VERSION VERSION_LESS 3.12)
      list(APPEND MDI_DEP_LIBS ${PYTHON_LIBRARIES})
    else()
      list(APPEND MDI_DEP_LIBS Python::Python)
    endif()

  endif()
  # need to add support for dlopen/dlsym, except when compiling for Windows.
  if(NOT (CMAKE_SYSTEM_NAME STREQUAL "Windows"))
    list(APPEND MDI_DEP_LIBS "${CMAKE_DL_LIBS}")
  endif()
  if(MDI_DEP_LIBS)
    set_target_properties(LAMMPS::MDI PROPERTIES
      IMPORTED_LINK_INTERFACE_LIBRARIES "${MDI_DEP_LIBS}")
  endif()

  target_link_libraries(lammps PRIVATE LAMMPS::MDI)
  target_link_libraries(lmp PRIVATE LAMMPS::MDI)

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
  target_include_directories(lmp PRIVATE ${mdi_INCLUDE_DIR})
  target_link_libraries(lmp PRIVATE ${mdi_LIBRARY})
endif()

target_compile_definitions(lammps PRIVATE -DLMP_MDI)
target_compile_definitions(lmp PRIVATE -DLMP_MDI)
