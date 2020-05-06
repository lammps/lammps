enable_language(C)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
  cmake_policy(SET CMP0074 NEW)
endif()
find_package(HDF5 REQUIRED)
target_link_libraries(h5md PRIVATE ${HDF5_LIBRARIES})
target_include_directories(h5md PUBLIC ${HDF5_INCLUDE_DIRS})
