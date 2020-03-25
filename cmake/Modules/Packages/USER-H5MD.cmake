if(PKG_USER-H5MD)
  enable_language(C)

  find_package(HDF5 REQUIRED)
  target_link_libraries(h5md PRIVATE ${HDF5_LIBRARIES})
  target_include_directories(h5md PUBLIC ${HDF5_INCLUDE_DIRS})
endif()
