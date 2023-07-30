cmake_minimum_required(VERSION 3.10)

project(lammps-gui VERSION 0.9 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(LAMMPS_GUI_USE_PLUGIN "Load LAMMPS library dynamically at runtime" OFF)

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
  get_filename_component(LAMMPS_DIR ${LAMMPS_SOURCE_DIR}/.. ABSOLUTE)
  # find LAMMPS library
  find_library(LAMMPS_LIBRARY
    NAMES lammps lammps_serial
    HINTS ${LAMMPS_LIB_DIR} ${LAMMPS_HEADER_DIR}/../build ${LAMMPS_HEADER_DIR}/../build-clang ${LAMMPS_HEADER_DIR}/../build-test ${LAMMPS_HEADER_DIR}
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
  
set(LAMMPS_PLUGINLIB_DIR ${LAMMPS_DIR}/examples/COUPLE/plugin)
find_package(QT NAMES Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)

if(LAMMPS_GUI_USE_PLUGIN)
  enable_language(C)
  set(PLUGIN_LOADER_SRC ${LAMMPS_PLUGINLIB_DIR}/liblammpsplugin.c)
endif()

set(PROJECT_SOURCES
  main.cpp
  codeeditor.cpp
  codeeditor.h
  highlighter.cpp
  highlighter.h
  linenumberarea.h
  lammpsgui.cpp
  lammpsgui.h
  lammpsgui.ui
  lammpsrunner.cpp
  lammpsrunner.h
  preferences.cpp
  preferences.h
  stdcapture.cpp
  ${PLUGIN_LOADER_SRC}
)
qt_add_resources(PROJECT_SOURCES lammpsgui.qrc)

if(APPLE)
  set(MACOSX_ICON_FILE ${LAMMPS_DIR}/cmake/packaging/lammps.icns)
  set(MACOSX_README_FILE ${LAMMPS_DIR}/cmake/packaging/README.macos)
  set(MACOSX_BACKGROUND_FILE ${LAMMPS_DIR}/cmake/packaging/LAMMPS_DMG_Background.png)
endif()

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(lammps-gui
        MANUAL_FINALIZATION
        ${MACOS_ICON_FILE}
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
            ${MACOSX_ICON_FILE}
            ${PROJECT_SOURCES}
        )
    endif()
endif()

if(LAMMPS_GUI_USE_PLUGIN)
  target_compile_definitions(lammps-gui PRIVATE LAMMPS_GUI_USE_PLUGIN)
  target_include_directories(lammps-gui PRIVATE ${LAMMPS_PLUGINLIB_DIR})
else()
  target_include_directories(lammps-gui PRIVATE ${LAMMPS_HEADER_DIR})
  target_link_libraries(lammps-gui PRIVATE ${LAMMPS_LIBRARY})
endif()
target_compile_definitions(lammps-gui PRIVATE LAMMPS_GUI_VERSION="${PROJECT_VERSION}")
target_link_libraries(lammps-gui PRIVATE Qt${QT_VERSION_MAJOR}::Widgets)
if(BUILD_OMP)
  find_package(OpenMP COMPONENTS CXX REQUIRED)
  target_link_libraries(lammps-gui PRIVATE OpenMP::OpenMP_CXX)
endif()

if(APPLE)
  set_target_properties(lammps-gui PROPERTIES
    MACOSX_BUNDLE_INFO_PLIST ${LAMMPS_DIR}/cmake/packaging/MacOSXBundleInfo.plist.in
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE_ICON_FILE lammps.icns
    MACOSX_BUNDLE_COPYRIGHT "(c) 2003 - 2023, The LAMMPS Developers"
    MACOSX_BUNDLE TRUE
  )
  # additional steps to populate the bundle tree and create the .dmg image file
  set(APP_CONTENTS ${CMAKE_BINARY_DIR}/lammps-gui.app/Contents)
  set_custom_target(complete-bundle ALL
    ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/bin
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/lmp ${APP_CONTENTS}/bin/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/msi2lmp ${APP_CONTENTS}/bin/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/binary2txt ${APP_CONTENTS}/bin/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/stl_bin2txt ${APP_CONTENTS}/bin/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/phana ${APP_CONTENTS}/bin/
    COMMAND ${CMAKE_COMMAND} -E create_symlink ../MacOS/lammps-gui ${APP_CONTENTS}/bin/lammps-gui
    COMMAND ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/Resources
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/cmake/packaging/README.macos ${APP_CONTENTS}/Resources/README.txt
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/cmake/packaging/LAMMPS_DMG_Background.png ${APP_CONTENTS}/Resources
    COMMAND ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/share/lammps
    COMMAND ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/share/lammps/man/man1
    COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${LAMMPS_DIR}/potentials ${APP_CONTENTS}/share/lammps/
    COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${LAMMPS_DIR}/bench ${APP_CONTENTS}/share/lammps/
    COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${LAMMPS_DIR}/tools/msi2lmp/frc_files ${APP_CONTENTS}/share/lammps/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/doc/lammps.1 ${APP_CONTENTS}/share/lammps/man/man1/
    COMMAND ${CMAKE_COMMAND} -E create_symlink lammps.1 ${APP_CONTENTS}/share/lammps/man/man1/lmp.1
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/doc/msi2lmp.1 ${APP_CONTENTS}/share/lammps/man/man1
    DEPENDS lammps-gui lammps
    COMMENT "Copying additional files into macOS app bundle tree"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  set_custom_target(build-dmg ALL
    COMMAND ${LAMMPS_DIR}/cmake/packaging/build-macos-dmg-with-gui.sh
    DEPENDS complete-bundle
    COMMENT "Create Drag-n-Drop installer disk image from app bundle"
    BYPRODUCT LAMMPS-macOS-multiarch.dmg
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(lammps-gui)
endif()
