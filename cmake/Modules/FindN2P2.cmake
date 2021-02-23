include(FindPackageHandleStandardArgs)

if (DEFINED ENV{N2P2_DIR})
  set(N2P2_DIR "${N2P2_DIR}")
endif()
message(STATUS "N2P2_DIR=${N2P2_DIR}")

find_path(N2P2_INCLUDE_DIR InterfaceLammps.h PATHS "${N2P2_DIR}/include")
find_library(N2P2_LIBNNP NAMES nnp PATHS "${N2P2_DIR}/lib")
find_library(N2P2_LIBNNPIF NAMES nnpif PATHS "${N2P2_DIR}/lib")

find_package_handle_standard_args(N2P2 DEFAULT_MSG N2P2_INCLUDE_DIR N2P2_LIBNNP)

if(N2P2_FOUND)
  set(N2P2_INCLUDE_DIRS ${N2P2_INCLUDE_DIR})
  set(N2P2_LIBRARIES ${N2P2_LIBNNPIF} ${N2P2_LIBNNP})

  mark_as_advanced(
    N2P2_DIR
    N2P2_INCLUDE_DIR
    N2P2_LIBNNP
    N2P2_LIBNNPIF
  )
endif()

