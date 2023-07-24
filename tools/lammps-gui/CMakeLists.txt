cmake_minimum_required(VERSION 3.10)

project(lammps-gui VERSION 0.1 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# when this file is included as subdirectory in the LAMMPS build, many settings are directly imported
if(LAMMPS_DIR)
  set(LAMMPS_HEADER_DIR ${LAMMPS_SOURCE_DIR})
  set(LAMMPS_LIBRARY lammps)
else()
  # NOTE: the next line should be commented out when used outside of the LAMMPS package
  get_filename_component(LAMMPS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/../../src ABSOLUTE)
  set(LAMMPS_HEADER_DIR ${LAMMPS_SOURCE_DIR} CACHE PATH "Location of LAMMPS headers")
  if(NOT LAMMPS_HEADER_DIR)
    message(FATAL_ERROR "Must set LAMMPS_HEADER_DIR")
  endif()
  # find LAMMPS library
  find_library(LAMMPS_LIBRARY
    NAMES lammps lammps_serial
    HINTS ${LAMMPS_HEADER_DIR}/../build ${LAMMPS_HEADER_DIR}/../build-clang ${LAMMPS_HEADER_DIR}/../build-test ${LAMMPS_HEADER_DIR}
    REQUIRED
  )
  # by default, install into $HOME/.local (not /usr/local),
  # so that no root access (and sudo) is needed
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Default install path" FORCE)
  endif()
  # ugly hacks for MSVC which by default always reports an old C++ standard in the __cplusplus macro
  # and prints lots of pointless warnings about "unsafe" functions
  if(MSVC)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
      add_compile_options(/Zc:__cplusplus)
      add_compile_options(/wd4244)
      add_compile_options(/wd4267)
      if(LAMMPS_EXCEPTIONS)
        add_compile_options(/EHsc)
      endif()
    endif()
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
  endif()
endif()

find_package(QT NAMES Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)

set(PROJECT_SOURCES
        main.cpp
        lammpsgui.cpp
        lammpsgui.h
        lammpsgui.ui
)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(lammps-gui
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
    )
# Define target properties for Android with Qt 6 as:
#    set_property(TARGET lammps-gui APPEND PROPERTY QT_ANDROID_PACKAGE_SOURCE_DIR
#                 ${CMAKE_CURRENT_SOURCE_DIR}/android)
# For more information, see https://doc.qt.io/qt-6/qt-add-executable.html#target-creation
else()
    if(ANDROID)
        add_library(lammps-gui SHARED
            ${PROJECT_SOURCES}
        )
# Define properties for Android with Qt 5 after find_package() calls as:
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
    else()
        add_executable(lammps-gui
            ${PROJECT_SOURCES}
        )
    endif()
endif()

target_include_directories(lammps-gui PRIVATE ${LAMMPS_HEADER_DIR})
target_link_libraries(lammps-gui PRIVATE Qt${QT_VERSION_MAJOR}::Widgets)
target_link_libraries(lammps-gui PRIVATE ${LAMMPS_LIBRARY})

set_target_properties(lammps-gui PROPERTIES
    MACOSX_BUNDLE_GUI_IDENTIFIER my.example.com
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(lammps-gui)
endif()
