enable_language(C)

# silence warning about ignoring <PackageName>_ROOT and use it.
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

find_package(HDF5 REQUIRED)
target_link_libraries(h5md PRIVATE ${HDF5_LIBRARIES})
target_include_directories(h5md PUBLIC ${HDF5_INCLUDE_DIRS})
