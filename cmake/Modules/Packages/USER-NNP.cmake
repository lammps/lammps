find_package(N2P2 REQUIRED)
target_include_directories(lammps PRIVATE ${N2P2_INCLUDE_DIRS})
target_link_libraries(lammps PRIVATE ${N2P2_LIBRARIES})
include(${N2P2_CMAKE_EXTRAS})
