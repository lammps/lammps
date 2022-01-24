include(FindPackageHandleStandardArgs)

# Check if N2P2_DIR is set manually.
if (DEFINED ENV{N2P2_DIR})
  set(N2P2_DIR "${N2P2_DIR}")
# If not, try if directory "lib/hdnnp/n2p2" exists.
else()
  get_filename_component(_fullpath "${LAMMPS_LIB_SOURCE_DIR}/hdnnp/n2p2" REALPATH)
  if (EXISTS ${_fullpath})
    set(N2P2_DIR "${_fullpath}")
  endif()
endif()

# Set path to include directory.
find_path(N2P2_INCLUDE_DIR InterfaceLammps.h HINTS "${N2P2_DIR}/include")
# Two libraries need to be linked: libnnp and libnnpif.
find_library(N2P2_LIBNNP NAMES nnp HINTS "${N2P2_DIR}/lib")
find_library(N2P2_LIBNNPIF NAMES nnpif HINTS "${N2P2_DIR}/lib")
# Users could compile n2p2 with extra flags which are then also required for
# pair_hdnnp.cpp compilation. To forward them to the LAMMPS build process n2p2
# writes a file with cmake commands, e.g.
#
# target_compile_definitions(lammps PRIVATE -DN2P2_NO_SF_GROUPS)
#
# to "lib/lammps-extra.cmake" which is then included by ML-HDNNP.cmake.
find_file(N2P2_CMAKE_EXTRA NAMES lammps-extra.cmake HINTS "${N2P2_DIR}/lib")

find_package_handle_standard_args(N2P2 DEFAULT_MSG
  N2P2_DIR
  N2P2_INCLUDE_DIR
  N2P2_LIBNNP
  N2P2_LIBNNPIF
  N2P2_CMAKE_EXTRA)

if(N2P2_FOUND)
  if (NOT TARGET N2P2::N2P2)
    # n2p2 core library "libnnp"
    add_library(N2P2::LIBNNP UNKNOWN IMPORTED)
    set_target_properties(N2P2::LIBNNP PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES ${N2P2_INCLUDE_DIR}
      IMPORTED_LOCATION ${N2P2_LIBNNP})
    # n2p2 interface library "libnnpif"
    add_library(N2P2::LIBNNPIF UNKNOWN IMPORTED)
    set_target_properties(N2P2::LIBNNPIF PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES ${N2P2_INCLUDE_DIR}
      IMPORTED_LOCATION ${N2P2_LIBNNPIF})
    # Put libnnp, libnnpif and include directory together.
    add_library(N2P2::N2P2 INTERFACE IMPORTED)
    set_property(TARGET N2P2::N2P2 PROPERTY
      INTERFACE_LINK_LIBRARIES N2P2::LIBNNPIF N2P2::LIBNNP)
    set(N2P2_CMAKE_EXTRAS ${N2P2_CMAKE_EXTRA})
  endif()
endif()

mark_as_advanced(
  N2P2_DIR
  N2P2_INCLUDE_DIR
  N2P2_LIBNNP
  N2P2_LIBNNPIF
  N2P2_CMAKE_EXTRA
)
