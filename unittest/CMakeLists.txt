
include(ExternalProject)
ExternalProject_Add(gtest
                    GIT_REPOSITORY  https://github.com/google/googletest.git
                    GIT_TAG         release-1.10.0
                    SOURCE_DIR      "${CMAKE_BINARY_DIR}/gtest-src"
                    BINARY_DIR      "${CMAKE_BINARY_DIR}/gtest-build"
                    CMAKE_ARGS      ${CMAKE_REQUEST_PIC} ${CMAKE_EXTRA_GTEST_OPTS}
                                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                                    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                                    -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
                                    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
                    BUILD_COMMAND   ${CMAKE_COMMAND} --build .
                    BUILD_BYPRODUCTS <BINARY_DIR>/lib/libgtest.a
                    INSTALL_COMMAND ""
                    TEST_COMMAND    "")

ExternalProject_Get_Property(gtest SOURCE_DIR)
set(GTEST_INCLUDE_DIR ${SOURCE_DIR}/googletest/include)
ExternalProject_Get_Property(gtest BINARY_DIR)
set(GTEST_LIBRARY_PATH ${BINARY_DIR}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

include(FindThreads)
set(GTEST_LIBRARY GTest::GTest)
add_library(${GTEST_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GTEST_LIBRARY} PROPERTIES
        IMPORTED_LOCATION ${GTEST_LIBRARY_PATH}
        INTERFACE_INCLUDE_DIRECTORIES ${GTEST_INCLUDE_DIR}
        INTERFACE_LINK_LIBRARIES Threads::Threads)
add_dependencies(${GTEST_LIBRARY} gtest)

add_subdirectory(cplusplus)
add_subdirectory(c-library)
add_subdirectory(fortran)
add_subdirectory(python)
