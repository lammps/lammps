enable_language(C)

add_library(qmmm STATIC ${LAMMPS_LIB_SOURCE_DIR}/qmmm/libqmmm.c)
set_target_properties(qmmm PROPERTIES OUTPUT_NAME lammps_qmmm${LAMMPS_MACHINE})
target_link_libraries(lammps PRIVATE qmmm)
target_include_directories(qmmm PUBLIC ${LAMMPS_LIB_SOURCE_DIR}/qmmm)
