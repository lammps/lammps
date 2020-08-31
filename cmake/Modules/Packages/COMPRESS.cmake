find_package(ZLIB REQUIRED)
target_link_libraries(lammps PRIVATE ZLIB::ZLIB)

find_package(Zstd)

if(Zstd_FOUND)
    target_compile_definitions(lammps PRIVATE -DLAMMPS_ZSTD)
    target_link_libraries(lammps PRIVATE Zstd::Zstd)
endif()
