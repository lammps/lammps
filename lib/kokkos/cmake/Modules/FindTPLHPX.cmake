
FIND_PACKAGE(HPX REQUIRED)
#as of right now, HPX doesn't export correctly
#so let's convert it to an interface target
KOKKOS_CREATE_IMPORTED_TPL(HPX INTERFACE
  LINK_LIBRARIES ${HPX_LIBRARIES}
  INCLUDES ${HPX_INCLUDE_DIRS}
)
#this is a bit funky since this is a CMake target
#but HPX doesn't export itself correctly
KOKKOS_EXPORT_CMAKE_TPL(HPX)

#I would prefer all of this gets replaced with
#KOKKOS_IMPORT_CMAKE_TPL(HPX)

