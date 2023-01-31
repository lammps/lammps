find_package(N2P2 QUIET)
if(N2P2_FOUND)
  set(DOWNLOAD_N2P2_DEFAULT OFF)
else()
  set(DOWNLOAD_N2P2_DEFAULT ON)
endif()
option(DOWNLOAD_N2P2 "Download n2p2 library instead of using an already installed one)" ${DOWNLOAD_N2P2_DEFAULT})
if(DOWNLOAD_N2P2)
  set(N2P2_URL "https://github.com/CompPhysVienna/n2p2/archive/v2.1.4.tar.gz" CACHE STRING "URL for n2p2 tarball")
  set(N2P2_MD5 "9595b066636cd6b90b0fef93398297a5" CACHE STRING "MD5 checksum of N2P2 tarball")
  mark_as_advanced(N2P2_URL)
  mark_as_advanced(N2P2_MD5)
  GetFallbackURL(N2P2_URL N2P2_FALLBACK)

  # adjust settings from detected compiler to compiler platform in n2p2 library
  # set compiler specific flag to turn on C++11 syntax (required on macOS and CentOS 7)
  if((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
    set(N2P2_COMP llvm)
    set(N2P2_CXX_STD "-std=c++11")
  elseif((CMAKE_CXX_COMPILER_ID STREQUAL "Intel") OR (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM"))
    set(N2P2_COMP intel)
    set(N2P2_CXX_STD "-std=c++11")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(N2P2_COMP gnu)
    set(N2P2_CXX_STD "-std=gnu++11")
  elseif((CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR (CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC"))
    set(N2P2_COMP gnu)
    set(N2P2_CXX_STD "--c++11")
  else() # default
    set(N2P2_COMP "")
  endif()

  # pass on archive creator command. prefer compiler specific version, if set.
  # important when using cross compiler.
  if(CMAKE_CXX_COMPILER_AR)
    set(N2P2_AR ${CMAKE_CXX_COMPILER_AR})
  else()
    set(N2P2_AR ${CMAKE_AR})
  endif()

  # adjust compilation of n2p2 library to whether MPI is requested in LAMMPS or not
  # need special care for compiling for MPICH2 with Linux-to-Windows cross compiler.
  if(NOT BUILD_MPI)
    set(N2P2_PROJECT_OPTIONS "-DN2P2_NO_MPI")
  else()
    # get path to MPI include directory
    get_target_property(N2P2_MPI_INCLUDE MPI::MPI_CXX INTERFACE_INCLUDE_DIRECTORIES)
    foreach (_INCL ${N2P2_MPI_INCLUDE})
      set(N2P2_PROJECT_OPTIONS "${N2P2_PROJECT_OPTIONS} -I${_INCL}")
    endforeach()
  endif()

  # prefer GNU make, if available. N2P2 lib seems to need it.
  find_program(N2P2_MAKE NAMES gmake make)

  # override compiler (optimization) flags in n2p2 library to flags used for LAMMPS
  # specifically -march=native can result in problems when compiling on HPC clusters or with a cross compiler
  # this convoluted way gets correct quoting/escaping when configuring the external project
  string(TOUPPER "${CMAKE_BUILD_TYPE}" BTYPE)
  set(N2P2_BUILD_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BTYPE}} ${N2P2_CXX_STD}")
  set(N2P2_BUILD_OPTIONS INTERFACES=LAMMPS COMP=${N2P2_COMP} "PROJECT_OPTIONS=${N2P2_PROJECT_OPTIONS}" "PROJECT_DEBUG="
    "PROJECT_CC=${CMAKE_CXX_COMPILER}" "PROJECT_MPICC=${CMAKE_CXX_COMPILER}" "PROJECT_CFLAGS=${N2P2_BUILD_FLAGS}"
    "PROJECT_AR=${N2P2_AR}" "APP_CORE=nnp-convert" "APP_TRAIN=nnp-train" "APP=nnp-convert")
  # echo final flag for debugging
  message(STATUS "N2P2 BUILD OPTIONS: ${N2P2_BUILD_OPTIONS}")

  # must have "sed" command to compile n2p2 library (for now)
  find_program(HAVE_SED sed)
  if(NOT HAVE_SED)
    message(FATAL_ERROR "Must have 'sed' program installed to compile 'n2p2' library for ML-HDNNP package")
  endif()

  # download compile n2p2 library. much patch MPI calls in LAMMPS interface to accommodate MPI-2 (e.g. for cross-compiling)
  include(ExternalProject)
  ExternalProject_Add(n2p2_build
    URL     ${N2P2_URL} ${N2P2_FALLBACK}
    URL_MD5 ${N2P2_MD5}
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    PATCH_COMMAND sed -i -e "s/\\(MPI_\\(P\\|Unp\\)ack(\\)/\\1(void *) /" src/libnnpif/LAMMPS/InterfaceLammps.cpp
    BUILD_COMMAND ${N2P2_MAKE} -C <SOURCE_DIR>/src -f makefile libnnpif ${N2P2_BUILD_OPTIONS}
    BUILD_ALWAYS YES
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    LOG_BUILD ON
    SOURCE_SUBDIR src/
    BUILD_BYPRODUCTS <SOURCE_DIR>/lib/libnnp.a <SOURCE_DIR>/lib/libnnpif.a
    )

  # create imported target LAMMPS::N2P2 from two libraries nnp and nnpif
  ExternalProject_get_property(n2p2_build SOURCE_DIR)
  # n2p2 core library "libnnp"
  add_library(LAMMPS::N2P2::LIBNNP UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::N2P2::LIBNNP PROPERTIES
    IMPORTED_LOCATION "${SOURCE_DIR}/lib/libnnp.a"
    INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/include")
  # n2p2 interface library "libnnpif"
  add_library(LAMMPS::N2P2::LIBNNPIF UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::N2P2::LIBNNPIF PROPERTIES
    IMPORTED_LOCATION "${SOURCE_DIR}/lib/libnnpif.a"
    INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/include")
  # nnpif library has MPI calls if MPI is enabled, so we must link with MPI libs
  if(BUILD_MPI)
    set_target_properties(LAMMPS::N2P2::LIBNNPIF PROPERTIES
      INTERFACE_LINK_LIBRARIES MPI::MPI_CXX)
    if((CMAKE_SYSTEM_NAME STREQUAL Windows) AND CMAKE_CROSSCOMPILING)
      add_dependencies(LAMMPS::N2P2::LIBNNPIF MPI::MPI_CXX)
    endif()
  endif()

  # final step to define imported target
  add_library(LAMMPS::N2P2 INTERFACE IMPORTED)
  set_property(TARGET LAMMPS::N2P2 PROPERTY
    INTERFACE_LINK_LIBRARIES LAMMPS::N2P2::LIBNNPIF LAMMPS::N2P2::LIBNNP)
  target_link_libraries(lammps PRIVATE LAMMPS::N2P2)

  add_dependencies(LAMMPS::N2P2 n2p2_build)
  # work around issues with older CMake versions
  file(MAKE_DIRECTORY "${SOURCE_DIR}/include")
  file(MAKE_DIRECTORY "${SOURCE_DIR}/lib")
else()
  find_package(N2P2)
  if(NOT N2P2_FOUND)
    message(FATAL_ERROR "n2p2 not found, help CMake to find it by setting N2P2_DIR, or set DOWNLOAD_N2P2=ON to download it")
  endif()
  target_link_libraries(lammps PRIVATE N2P2::N2P2)
  include(${N2P2_CMAKE_EXTRAS})
endif()
