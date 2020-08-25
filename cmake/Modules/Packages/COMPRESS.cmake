find_package(ZLIB REQUIRED)
target_link_libraries(lammps PRIVATE ZLIB::ZLIB)

find_package(Zstd REQUIRED)
target_link_libraries(lammps PRIVATE Zstd::Zstd)
