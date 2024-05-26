cmake_minimum_required(VERSION 3.16)

project(lammps-gui VERSION 1.5.11 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(LAMMPS_GUI_USE_PLUGIN "Load LAMMPS library dynamically at runtime" OFF)
mark_as_advanced(LAMMPS_GUI_USE_PLUGIN)
option(LAMMPS_GUI_USE_QT5 "Prefer using Qt5 over Qt6" OFF)

include(CheckIncludeFileCXX)
# helper function to check for usable omp.h header
function(check_omp_h_include)
  find_package(OpenMP COMPONENTS CXX QUIET)
  if(OpenMP_CXX_FOUND)
    set(CMAKE_REQUIRED_FLAGS ${OpenMP_CXX_FLAGS})
    set(CMAKE_REQUIRED_INCLUDES ${OpenMP_CXX_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LINK_OPTIONS ${OpenMP_CXX_FLAGS})
    set(CMAKE_REQUIRED_LIBRARIES ${OpenMP_CXX_LIBRARIES})
    check_include_file_cxx(omp.h _have_omp_h)
  else()
    set(_have_omp_h FALSE)
  endif()
  set(HAVE_OMP_H_INCLUDE ${_have_omp_h} PARENT_SCOPE)
endfunction()

# detect if we may enable OpenMP support by default
set(BUILD_OMP_DEFAULT OFF)
find_package(OpenMP COMPONENTS CXX QUIET)
if(OpenMP_CXX_FOUND)
  check_omp_h_include()
  if(HAVE_OMP_H_INCLUDE)
    set(BUILD_OMP_DEFAULT ON)
  endif()
endif()

option(BUILD_OMP "Build with OpenMP support" ${BUILD_OMP_DEFAULT})
if(BUILD_OMP)
  find_package(OpenMP COMPONENTS CXX REQUIRED)
  check_omp_h_include()
  if(NOT HAVE_OMP_H_INCLUDE)
    message(FATAL_ERROR "Cannot find the 'omp.h' header file required for full OpenMP support")
  endif()
endif()

# checks
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
      add_compile_options(/wd4250)
      if(LAMMPS_EXCEPTIONS)
        add_compile_options(/EHsc)
      endif()
    endif()
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
  endif()
endif()

#get_lammps_version(${LAMMPS_SOURCE_DIR}/version.h LAMMPS_VERSION)
#if (LAMMPS_VERSION LESS_EQUAL 20230802)
#  message(FATAL_ERROR "LAMMPS GUI v${PROJECT_VERSION} is incompatible with LAMMPS 2 August 2023 or older")
#endif()

set(LAMMPS_PLUGINLIB_DIR ${LAMMPS_DIR}/examples/COUPLE/plugin)
if(LAMMPS_GUI_USE_PLUGIN)
  if(APPLE)
    message(FATAL_ERROR "Dynamic loading of the LAMMPS library currently not supported on this platform")
  endif()
  enable_language(C)
  set(PLUGIN_LOADER_SRC ${LAMMPS_PLUGINLIB_DIR}/liblammpsplugin.c)
endif()

# include resource compiler to embed icons into the executable on Windows
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  enable_language(RC)
  set(ICON_RC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/lmpicons.rc)
endif()

# we require Qt 5 and at least version 5.12 at that.
if(NOT LAMMPS_GUI_USE_QT5)
  find_package(Qt6 6.2 QUIET COMPONENTS Widgets Charts)
endif()
if(NOT Qt6_FOUND)
  find_package(Qt5 5.12 REQUIRED COMPONENTS Widgets Charts)
  set(QT_VERSION_MAJOR 5)
else()
  set(QT_VERSION_MAJOR 6)
endif()
message(STATUS "Using Qt version ${Qt${QT_VERSION_MAJOR}_VERSION} for LAMMPS GUI")

set(PROJECT_SOURCES
  main.cpp
  lammpsgui.cpp
  lammpsgui.h
  lammpsgui.ui
  chartviewer.cpp
  chartviewer.h
  codeeditor.cpp
  codeeditor.h
  helpers.cpp
  highlighter.cpp
  highlighter.h
  imageviewer.cpp
  imageviewer.h
  lammpsrunner.h
  lammpswrapper.cpp
  lammpswrapper.h
  linenumberarea.h
  logwindow.cpp
  logwindow.h
  preferences.cpp
  preferences.h
  setvariables.cpp
  setvariables.h
  slideshow.h
  slideshow.cpp
  stdcapture.cpp
  ${PLUGIN_LOADER_SRC}
  ${ICON_RC_FILE}
)
if(QT_VERSION_MAJOR EQUAL 6)
  qt6_add_resources(PROJECT_SOURCES lammpsgui.qrc)
else()
  qt5_add_resources(PROJECT_SOURCES lammpsgui.qrc)
endif()

if(APPLE)
  set(MACOSX_ICON_FILE ${LAMMPS_DIR}/cmake/packaging/lammps.icns)
  set(MACOSX_README_FILE ${LAMMPS_DIR}/cmake/packaging/README.macos)
  set(MACOSX_BACKGROUND_FILE ${LAMMPS_DIR}/cmake/packaging/LAMMPS_DMG_Background.png)
endif()

if(QT_VERSION_MAJOR EQUAL 6)
  qt_add_executable(lammps-gui
    MANUAL_FINALIZATION
    ${MACOSX_ICON_FILE}
    ${PROJECT_SOURCES}
  )
else()
  add_executable(lammps-gui
    ${MACOSX_ICON_FILE}
    ${PROJECT_SOURCES}
  )
endif()

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(lammps-gui)
endif()

# compilation settings
if(LAMMPS_GUI_USE_PLUGIN)
  target_compile_definitions(lammps-gui PRIVATE LAMMPS_GUI_USE_PLUGIN)
  target_include_directories(lammps-gui PRIVATE ${LAMMPS_PLUGINLIB_DIR})
else()
  target_include_directories(lammps-gui PRIVATE ${LAMMPS_HEADER_DIR})
  target_link_libraries(lammps-gui PRIVATE ${LAMMPS_LIBRARY})
endif()
target_include_directories(lammps-gui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(lammps-gui PRIVATE LAMMPS_GUI_VERSION="${PROJECT_VERSION}")
target_link_libraries(lammps-gui PRIVATE Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::Charts)
if(BUILD_OMP)
  find_package(OpenMP COMPONENTS CXX REQUIRED)
  target_link_libraries(lammps-gui PRIVATE OpenMP::OpenMP_CXX)
endif()

# when compiling on macOS, create an "app bundle"
if(APPLE)
  set_target_properties(lammps-gui PROPERTIES
    MACOSX_BUNDLE_INFO_PLIST ${LAMMPS_DIR}/cmake/packaging/MacOSXBundleInfo.plist.in
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE_ICON_FILE lammps.icns
    MACOSX_BUNDLE_COPYRIGHT "(c) 2003 - 2023, The LAMMPS Developers"
    MACOSX_BUNDLE TRUE
  )
  # additional targets to populate the bundle tree and create the .dmg image file
  set(APP_CONTENTS ${CMAKE_BINARY_DIR}/lammps-gui.app/Contents)
  add_custom_target(complete-bundle
    ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/bin
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/lmp ${APP_CONTENTS}/bin/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/msi2lmp ${APP_CONTENTS}/bin/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/binary2txt ${APP_CONTENTS}/bin/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/stl_bin2txt ${APP_CONTENTS}/bin/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/phana ${APP_CONTENTS}/bin/
    COMMAND ${CMAKE_COMMAND} -E create_symlink ../MacOS/lammps-gui ${APP_CONTENTS}/bin/lammps-gui
    COMMAND ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/Resources
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/cmake/packaging/README.macos ${APP_CONTENTS}/Resources/README.txt
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/cmake/packaging/lammps.icns ${APP_CONTENTS}/Resources
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/cmake/packaging/LAMMPS_DMG_Background.png ${APP_CONTENTS}/Resources
    COMMAND ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/share/lammps
    COMMAND ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/share/lammps/man/man1
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${LAMMPS_DIR}/potentials ${APP_CONTENTS}/share/lammps/potentials
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${LAMMPS_DIR}/bench ${APP_CONTENTS}/share/lammps/bench
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${LAMMPS_DIR}/tools/msi2lmp/frc_files ${APP_CONTENTS}/share/lammps/frc_files
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/doc/lammps.1 ${APP_CONTENTS}/share/lammps/man/man1/
    COMMAND ${CMAKE_COMMAND} -E create_symlink lammps.1 ${APP_CONTENTS}/share/lammps/man/man1/lmp.1
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/doc/msi2lmp.1 ${APP_CONTENTS}/share/lammps/man/man1
    DEPENDS lammps-gui lammps lmp binary2txt stl_bin2txt msi2lmp phana
    COMMENT "Copying additional files into macOS app bundle tree"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  if(FFMPEG_EXECUTABLE)
    add_custom_target(copy-ffmpeg
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${FFMPEG_EXECUTABLE} ${APP_CONTENTS}/bin/
      COMMENT "Copying FFMpeg into macOS app bundle tree"
      DEPENDS complete-bundle
    )
    set(FFMPEG_TARGET copy-ffmpeg)
  endif()
  add_custom_target(dmg
    COMMAND ${LAMMPS_DIR}/cmake/packaging/build_macos_dmg.sh
    DEPENDS complete-bundle ${FFMPEG_TARGET}
    COMMENT "Create Drag-n-Drop installer disk image from app bundle"
    BYPRODUCT LAMMPS-macOS-multiarch.dmg
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  # settings or building on Windows with Visual Studio
elseif(MSVC)
  install(TARGETS lammps-gui DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES $<TARGET_RUNTIME_DLLS:lammps-gui> TYPE BIN)
  if(BUILD_SHARED_LIBS)
    install(FILES $<TARGET_RUNTIME_DLLS:lammps> TYPE BIN)
  endif()
  install(FILES $<TARGET_RUNTIME_DLLS:lmp> TYPE BIN)
  # find path to VC++ init batch file
  get_filename_component(VC_COMPILER_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)
  get_filename_component(VC_BASE_DIR "${VC_COMPILER_DIR}/../../../../../.." ABSOLUTE)
  set(VC_INIT "${VC_BASE_DIR}/Auxiliary/Build/vcvarsall.bat")
  get_filename_component(QT5_BIN_DIR "${Qt5Core_DIR}/../../../bin" ABSOLUTE)
  get_filename_component(INSTNAME ${CMAKE_INSTALL_PREFIX} NAME)
  install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -D INSTNAME=${INSTNAME} -D VC_INIT=\"${VC_INIT}\" -D QT5_BIN_DIR=\"${QT5_BIN_DIR}\" -P \"${CMAKE_SOURCE_DIR}/packaging/build_windows_vs.cmake\" WORKING_DIRECTORY \"${CMAKE_INSTALL_PREFIX}/..\" COMMAND_ECHO STDOUT)")
elseif((CMAKE_SYSTEM_NAME STREQUAL "Windows") AND CMAKE_CROSSCOMPILING)
  install(TARGETS lammps-gui DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES $<TARGET_RUNTIME_DLLS:lammps-gui> TYPE BIN)
  if(BUILD_SHARED_LIBS)
    install(FILES $<TARGET_RUNTIME_DLLS:lammps> TYPE BIN)
  endif()
  install(FILES $<TARGET_RUNTIME_DLLS:lmp> TYPE BIN)
  add_custom_target(zip
    COMMAND sh -vx ${LAMMPS_DIR}/cmake/packaging/build_windows_cross_zip.sh ${CMAKE_INSTALL_PREFIX}
    DEPENDS lmp lammps-gui
    COMMENT "Create zip file with windows binaries"
    BYPRODUCT LAMMPS-Win10-amd64.zip
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
elseif((CMAKE_SYSTEM_NAME STREQUAL "Linux") AND NOT LAMMPS_GUI_USE_PLUGIN)
  install(TARGETS lammps-gui DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/lammps-gui.desktop DESTINATION ${CMAKE_INSTALL_DATADIR}/applications/)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/lammps-input.xml DESTINATION ${CMAKE_INSTALL_DATADIR}/mime/packages/)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/lammps-input.xml DESTINATION ${CMAKE_INSTALL_DATADIR}/mime/text/x-application-lammps.xml)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/icons/hicolor DESTINATION ${CMAKE_INSTALL_DATADIR}/icons/)
  install(CODE [[
    file(GET_RUNTIME_DEPENDENCIES
      LIBRARIES $<TARGET_FILE:lammps>
      EXECUTABLES $<TARGET_FILE:lammps-gui> $<TARGET_FILE:lmp>
      RESOLVED_DEPENDENCIES_VAR _r_deps
      UNRESOLVED_DEPENDENCIES_VAR _u_deps
    )
    foreach(_file ${_r_deps})
      file(INSTALL
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        TYPE SHARED_LIBRARY
        FOLLOW_SYMLINK_CHAIN
        FILES "${_file}"
      )
    endforeach()
    list(LENGTH _u_deps _u_length)
    if("${_u_length}" GREATER 0)
      message(WARNING "Unresolved dependencies detected: ${_u_deps}")
    endif() ]]
  )
  add_custom_target(tgz
    COMMAND ${LAMMPS_DIR}/cmake/packaging/build_linux_tgz.sh
    DEPENDS lammps-gui
    COMMENT "Create compressed tar file of LAMMPS with dependent libraries and wrapper"
    BYPRODUCT LAMMPS-Linux-amd64.tar.gz
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()
