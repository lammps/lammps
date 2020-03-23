if(PKG_COMPRESS)
  find_package(ZLIB REQUIRED)
  target_link_libraries(lammps PRIVATE ZLIB::ZLIB)
endif()
