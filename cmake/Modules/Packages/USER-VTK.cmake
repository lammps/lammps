find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})
target_compile_definitions(lammps PRIVATE -DLAMMPS_VTK)
target_link_libraries(lammps PRIVATE ${VTK_LIBRARIES})
