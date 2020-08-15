include(GTest)

add_subdirectory(utils)
add_subdirectory(formats)
add_subdirectory(commands)
add_subdirectory(force-styles)

find_package(ClangFormat 8.0)

if(ClangFormat_FOUND)
  set(UNITTEST_SOURCES)
  file(GLOB_RECURSE UNITTEST_SOURCES *.cpp *.h)
  add_custom_target(format-tests
    COMMAND ${ClangFormat_EXECUTABLE} --verbose -i -style=file ${UNITTEST_SOURCES}
    DEPENDS ${UNITTEST_SOURCES})
endif()
