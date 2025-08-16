find_package(HPX REQUIRED 1.8.0)
#as of right now, HPX doesn't export correctly
#so let's convert it to an interface target
kokkos_create_imported_tpl(HPX INTERFACE LINK_LIBRARIES ${HPX_LIBRARIES} INCLUDES ${HPX_INCLUDE_DIRS})
#this is a bit funky since this is a CMake target
#but HPX doesn't export itself correctly
kokkos_export_cmake_tpl(HPX)

#I would prefer all of this gets replaced with
#KOKKOS_IMPORT_CMAKE_TPL(HPX)
