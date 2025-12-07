# -*- CMake -*- file for tests using DEM models

# build support library
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  string(REPLACE "/" "\\\\" TEST_INPUT_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/tests")
else()
  set(TEST_INPUT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endif()
add_library(granular_tests STATIC
  "${CMAKE_CURRENT_SOURCE_DIR}/../force-styles/error_stats.cpp"
  test_main.cpp)
#  test_config_reader.cpp
#  yaml_writer.cpp
if(YAML_FOUND)
  target_compile_definitions(granular_tests PRIVATE TEST_INPUT_FOLDER=${TEST_INPUT_FOLDER})
else()
  # we always use static linkage with local compiled libyaml
  target_compile_definitions(granular_tests PRIVATE TEST_INPUT_FOLDER=${TEST_INPUT_FOLDER} YAML_DECLARE_STATIC)
endif()
target_include_directories(granular_tests PRIVATE ${LAMMPS_SOURCE_DIR})
target_link_libraries(granular_tests PUBLIC GTest::GMock Yaml::Yaml lammps)

