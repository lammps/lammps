message(STATUS "Downloading and building Google Test library")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(GTEST_LIB_POSTFIX d)
else()
  set(GTEST_LIB_POSTFIX)
endif()

include(ExternalProject)
set(GTEST_URL "https://github.com/google/googletest/archive/release-1.10.0.tar.gz" CACHE STRING "URL for GTest tarball")
set(GTEST_MD5 "ecd1fa65e7de707cd5c00bdac56022cd" CACHE STRING "MD5 checksum of GTest tarball")
mark_as_advanced(GTEST_URL)
mark_as_advanced(GTEST_MD5)
ExternalProject_Add(googletest
                    URL             ${GTEST_URL}
                    URL_MD5         ${GTEST_MD5}
                    SOURCE_DIR      "${CMAKE_BINARY_DIR}/gtest-src"
                    BINARY_DIR      "${CMAKE_BINARY_DIR}/gtest-build"
                    CMAKE_ARGS      ${CMAKE_REQUEST_PIC} ${CMAKE_EXTRA_GTEST_OPTS}
                                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                                    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                                    -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
                                    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
                    BUILD_BYPRODUCTS <BINARY_DIR>/lib/libgtest${GTEST_LIB_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}
                                     <BINARY_DIR>/lib/libgmock${GTEST_LIB_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}
                                     <BINARY_DIR>/lib/libgtest_main${GTEST_LIB_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}
                                     <BINARY_DIR>/lib/libgmock_main${GTEST_LIB_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}
                    LOG_DOWNLOAD ON
                    LOG_CONFIGURE ON
                    LOG_BUILD ON
                    INSTALL_COMMAND ""
                    TEST_COMMAND    "")

ExternalProject_Get_Property(googletest SOURCE_DIR)
set(GTEST_INCLUDE_DIR ${SOURCE_DIR}/googletest/include)
set(GMOCK_INCLUDE_DIR ${SOURCE_DIR}/googlemock/include)

# workaround for CMake 3.10 on ubuntu 18.04
file(MAKE_DIRECTORY ${GTEST_INCLUDE_DIR})
file(MAKE_DIRECTORY ${GMOCK_INCLUDE_DIR})

ExternalProject_Get_Property(googletest BINARY_DIR)
set(GTEST_LIBRARY_PATH ${BINARY_DIR}/lib/libgtest${GTEST_LIB_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX})
set(GMOCK_LIBRARY_PATH ${BINARY_DIR}/lib/libgmock${GTEST_LIB_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX})
set(GTEST_MAIN_LIBRARY_PATH ${BINARY_DIR}/lib/libgtest_main${GTEST_LIB_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX})
set(GMOCK_MAIN_LIBRARY_PATH ${BINARY_DIR}/lib/libgmock_main${GTEST_LIB_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX})

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

find_package(Threads QUIET)

add_library(GTest::GTest UNKNOWN IMPORTED)
set_target_properties(GTest::GTest PROPERTIES
        IMPORTED_LOCATION ${GTEST_LIBRARY_PATH}
        INTERFACE_INCLUDE_DIRECTORIES ${GTEST_INCLUDE_DIR}
        INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
add_dependencies(GTest::GTest googletest)

add_library(GTest::GMock UNKNOWN IMPORTED)
set_target_properties(GTest::GMock PROPERTIES
        IMPORTED_LOCATION ${GMOCK_LIBRARY_PATH}
        INTERFACE_INCLUDE_DIRECTORIES ${GMOCK_INCLUDE_DIR}
        INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
add_dependencies(GTest::GMock googletest)

add_library(GTest::GTestMain UNKNOWN IMPORTED)
set_target_properties(GTest::GTestMain PROPERTIES
        IMPORTED_LOCATION ${GTEST_MAIN_LIBRARY_PATH}
        INTERFACE_INCLUDE_DIRECTORIES ${GTEST_INCLUDE_DIR}
        INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
add_dependencies(GTest::GTestMain googletest)

add_library(GTest::GMockMain UNKNOWN IMPORTED)
set_target_properties(GTest::GMockMain PROPERTIES
        IMPORTED_LOCATION ${GMOCK_MAIN_LIBRARY_PATH}
        INTERFACE_INCLUDE_DIRECTORIES ${GMOCK_INCLUDE_DIR}
        INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
add_dependencies(GTest::GMockMain googletest)
