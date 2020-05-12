message(STATUS "Downloading and building Google Test library")
include(ExternalProject)
ExternalProject_Add(googletest
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
                    BUILD_BYPRODUCTS <BINARY_DIR>/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a
                                     <BINARY_DIR>/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gmock.a
                                     <BINARY_DIR>/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a
                                     <BINARY_DIR>/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gmock_main.a
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
set(GTEST_LIBRARY_PATH ${BINARY_DIR}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)
set(GMOCK_LIBRARY_PATH ${BINARY_DIR}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gmock.a)
set(GTEST_MAIN_LIBRARY_PATH ${BINARY_DIR}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a)
set(GMOCK_MAIN_LIBRARY_PATH ${BINARY_DIR}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gmock_main.a)

set(GTEST_LIBRARY GTest::GTest)
set(GMOCK_LIBRARY GTest::GMock)
set(GTEST_MAIN_LIBRARY GTest::GTestMain)
set(GMOCK_MAIN_LIBRARY GTest::GMockMain)

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

find_package(Threads QUIET)

add_library(${GTEST_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GTEST_LIBRARY} PROPERTIES
        IMPORTED_LOCATION ${GTEST_LIBRARY_PATH}
        INTERFACE_INCLUDE_DIRECTORIES ${GTEST_INCLUDE_DIR}
        IMPORTED_LINK_INTERFACE_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
add_dependencies(${GTEST_LIBRARY} googletest)

add_library(${GMOCK_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GMOCK_LIBRARY} PROPERTIES
        IMPORTED_LOCATION ${GMOCK_LIBRARY_PATH}
        INTERFACE_INCLUDE_DIRECTORIES ${GMOCK_INCLUDE_DIR}
        IMPORTED_LINK_INTERFACE_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
add_dependencies(${GMOCK_LIBRARY} googletest)

add_library(${GTEST_MAIN_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GTEST_MAIN_LIBRARY} PROPERTIES
        IMPORTED_LOCATION ${GTEST_MAIN_LIBRARY_PATH}
        INTERFACE_INCLUDE_DIRECTORIES ${GTEST_INCLUDE_DIR}
        IMPORTED_LINK_INTERFACE_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
add_dependencies(${GTEST_MAIN_LIBRARY} googletest)

add_library(${GMOCK_MAIN_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GMOCK_MAIN_LIBRARY} PROPERTIES
        IMPORTED_LOCATION ${GMOCK_MAIN_LIBRARY_PATH}
        INTERFACE_INCLUDE_DIRECTORIES ${GMOCK_INCLUDE_DIR}
        IMPORTED_LINK_INTERFACE_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
add_dependencies(${GMOCK_MAIN_LIBRARY} googletest)
