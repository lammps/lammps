
add_executable(test_library_open test_library_open.cpp test_main.cpp)
target_link_libraries(test_library_open PRIVATE lammps GTest::GTest GTest::GMock)
add_test(LibraryOpen test_library_open)

add_executable(test_library_commands test_library_commands.cpp test_main.cpp)
target_link_libraries(test_library_commands PRIVATE lammps GTest::GTest GTest::GMock)
add_test(LibraryCommands test_library_commands)

add_executable(test_library_properties test_library_properties.cpp test_main.cpp)
target_link_libraries(test_library_properties PRIVATE lammps GTest::GTest GTest::GMock)
target_compile_definitions(test_library_properties PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
add_test(LibraryProperties test_library_properties)

set(PKG_COUNT 0)
foreach(PKG ${STANDARD_PACKAGES} ${SUFFIX_PACKAGES})
  if(PKG_${PKG})
    MATH(EXPR PKG_COUNT "${PKG_COUNT}+1")
  endif()
endforeach()

if(PKG_MANYBODY)
  set(HAS_MANYBODY 1)
else()
  set(HAS_MANYBODY 0)
endif()

add_executable(test_library_config test_library_config.cpp test_main.cpp)
target_link_libraries(test_library_config PRIVATE lammps GTest::GTest GTest::GMock)
target_compile_definitions(test_library_config PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR}
        -DNUM_LAMMPS_PACKAGES=${PKG_COUNT} -DLAMMPS_HAS_MANYBODY=${HAS_MANYBODY})
add_test(LibraryConfig test_library_config)
