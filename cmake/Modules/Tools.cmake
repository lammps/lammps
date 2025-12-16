if(BUILD_TOOLS)
  add_executable(binary2txt ${LAMMPS_TOOLS_DIR}/binary2txt.cpp)
  target_compile_definitions(binary2txt PRIVATE -DLAMMPS_${LAMMPS_SIZES})
  install(TARGETS binary2txt DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(stl_bin2txt ${LAMMPS_TOOLS_DIR}/stl_bin2txt.cpp)
  install(TARGETS stl_bin2txt DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(reformat-json ${LAMMPS_TOOLS_DIR}/json/reformat-json.cpp)
  target_include_directories(reformat-json PRIVATE ${LAMMPS_SOURCE_DIR})
  install(TARGETS reformat-json DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_custom_target(tools ALL COMMENT "Building tools")
  add_dependencies(tools binary2txt stl_bin2txt reformat-json)

  include(CheckGeneratorSupport)
  if(CMAKE_GENERATOR_SUPPORT_FORTRAN)
    include(CheckLanguage)
    check_language(Fortran)
    if(CMAKE_Fortran_COMPILER)
      enable_language(Fortran)
      add_executable(chain.x ${LAMMPS_TOOLS_DIR}/chain.f90)
      target_link_libraries(chain.x PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
      add_executable(micelle2d.x ${LAMMPS_TOOLS_DIR}/micelle2d.f90)
      target_link_libraries(micelle2d.x PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
      install(TARGETS chain.x micelle2d.x DESTINATION ${CMAKE_INSTALL_BINDIR})
      add_dependencies(tools chain.x micelle2d.x)
    else()
      message(WARNING "No suitable Fortran compiler found, skipping build of 'chain.x' and 'micelle2d.x'")
    endif()
  else()
    message(WARNING "CMake build doesn't support Fortran, skipping build of 'chain.x' and 'micelle2d.x'")
  endif()

  enable_language(C)
  get_filename_component(MSI2LMP_SOURCE_DIR ${LAMMPS_TOOLS_DIR}/msi2lmp/src ABSOLUTE)
  file(GLOB MSI2LMP_SOURCES CONFIGURE_DEPENDS ${MSI2LMP_SOURCE_DIR}/[^.]*.c)
  add_executable(msi2lmp ${MSI2LMP_SOURCES})
  if(STANDARD_MATH_LIB)
    target_link_libraries(msi2lmp PRIVATE ${STANDARD_MATH_LIB})
  endif()
  install(TARGETS msi2lmp DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES ${LAMMPS_DOC_DIR}/msi2lmp.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

  add_subdirectory(${LAMMPS_TOOLS_DIR}/phonon ${CMAKE_BINARY_DIR}/phana_build)
  add_dependencies(tools msi2lmp phana)
endif()

if(BUILD_LAMMPS_GUI)
  include(ExternalProject)
  # When building LAMMPS-GUI with LAMMPS we don't support plugin mode and don't include docs.
  ExternalProject_Add(lammps-gui_build
    GIT_REPOSITORY https://github.com/akohlmey/lammps-gui.git
    GIT_TAG main
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    CMAKE_ARGS -D BUILD_DOC=OFF
               -D LAMMPS_GUI_USE_PLUGIN=OFF
               -D LAMMPS_SOURCE_DIR=${LAMMPS_SOURCE_DIR}
               -D LAMMPS_LIBRARY=$<TARGET_FILE:lammps>
               -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>
               -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -D CMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
               -D CMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    DEPENDS lammps
    BUILD_BYPRODUCTS <INSTALL_DIR>/bin/lammps-gui
  )
  add_custom_target(lammps-gui ALL
          ${CMAKE_COMMAND} -E copy_if_different lammps-gui_build-prefix/bin/lammps-gui* ${CMAKE_BINARY_DIR}
          DEPENDS lammps-gui_build
  )

  # packaging support for LAMMPS-GUI when compiled with LAMMPS
  option(BUILD_WHAM "Download and compile WHAM executable from Grossfield Lab" YES)
  if(BUILD_WHAM)
    set(WHAM_URL "http://membrane.urmc.rochester.edu/sites/default/files/wham/wham-release-2.1.0.tgz" CACHE STRING "URL for WHAM tarball")
    set(WHAM_MD5 "4ed6e24254925ec124f44bb381c3b87f" CACHE STRING "MD5 checksum of WHAM tarball")
    mark_as_advanced(WHAM_URL)
    mark_as_advanced(WHAM_MD5)

    get_filename_component(archive ${WHAM_URL} NAME)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/_deps/src)
    if(EXISTS ${CMAKE_BINARY_DIR}/_deps/${archive})
      file(MD5 ${CMAKE_BINARY_DIR}/_deps/${archive} DL_MD5)
    endif()
    if(NOT "${DL_MD5}" STREQUAL "${WHAM_MD5}")
      message(STATUS "Downloading ${WHAM_URL}")
      file(DOWNLOAD ${WHAM_URL} ${CMAKE_BINARY_DIR}/_deps/${archive} STATUS DL_STATUS SHOW_PROGRESS)
      file(MD5 ${CMAKE_BINARY_DIR}/_deps/${archive} DL_MD5)
      if((NOT DL_STATUS EQUAL 0) OR (NOT "${DL_MD5}" STREQUAL "${WHAM_MD5}"))
        message(ERROR "Download of WHAM sources from ${WHAM_URL} failed")
      endif()
    else()
      message(STATUS "Using already downloaded archive ${CMAKE_BINARY_DIR}/_deps/${archive}")
    endif()
    message(STATUS "Unpacking and configuring ${archive}")

    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_BINARY_DIR}/_deps/${archive}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/_deps/src)
    find_package(Patch)
    if(PATCH_FOUND)
      message(STATUS "Apply patch to customize WHAM using ${Patch_EXECUTABLE}")
      execute_process(
        COMMAND ${Patch_EXECUTABLE} -p1 -i ${CMAKE_SOURCE_DIR}/cmake/packaging/update-wham-2.1.0.patch
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/_deps/src/wham/
      )
    endif()
    file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/_deps/wham-src)
    file(RENAME "${CMAKE_BINARY_DIR}/_deps/src/wham" ${CMAKE_BINARY_DIR}/_deps/wham-src)
    file(COPY packaging/CMakeLists.wham DESTINATION ${CMAKE_BINARY_DIR}/_deps/wham-src/)
    file(RENAME "${CMAKE_BINARY_DIR}/_deps/wham-src/CMakeLists.wham"
      "${CMAKE_BINARY_DIR}/_deps/wham-src/CMakeLists.txt")
    add_subdirectory("${CMAKE_BINARY_DIR}/_deps/wham-src" "${CMAKE_BINARY_DIR}/_deps/wham-build")
    set(WHAM_EXE wham wham-2d)
  endif()

  # build LAMMPS-GUI and LAMMPS as flatpak, if tools are installed
  find_program(FLATPAK_COMMAND flatpak DOC "Path to flatpak command")
  find_program(FLATPAK_BUILDER flatpak-builder DOC "Path to flatpak-builder command")
  if(FLATPAK_COMMAND AND FLATPAK_BUILDER)
    file(STRINGS ${LAMMPS_DIR}/src/version.h line REGEX LAMMPS_VERSION)
    string(REGEX REPLACE "#define LAMMPS_VERSION \"([0-9]+) ([A-Za-z][A-Za-z][A-Za-z])[A-Za-z]* ([0-9]+)\""
                        "\\1\\2\\3" LAMMPS_RELEASE "${line}")
    set(FLATPAK_BUNDLE "LAMMPS-Linux-x86_64-GUI-${LAMMPS_RELEASE}.flatpak")
    add_custom_target(flatpak
      COMMAND ${FLATPAK_COMMAND} --user remote-add --if-not-exists flathub https://dl.flathub.org/repo/flathub.flatpakrepo
      COMMAND ${FLATPAK_BUILDER} --force-clean --verbose --repo=${CMAKE_CURRENT_BINARY_DIR}/flatpak-repo
                               --install-deps-from=flathub --state-dir=${CMAKE_CURRENT_BINARY_DIR}
                               --user --ccache --default-branch=${LAMMPS_RELEASE}
                               flatpak-build ${LAMMPS_DIR}/cmake/packaging/org.lammps.lammps-gui.yml
      COMMAND ${FLATPAK_COMMAND} build-bundle --runtime-repo=https://flathub.org/repo/flathub.flatpakrepo --verbose
                               ${CMAKE_CURRENT_BINARY_DIR}/flatpak-repo
                               ${FLATPAK_BUNDLE} org.lammps.lammps-gui ${LAMMPS_RELEASE}
      COMMENT "Create Flatpak bundle file of LAMMPS and LAMMPS-GUI"
      BYPRODUCT ${FLATPAK_BUNDLE}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
  else()
    add_custom_target(flatpak
      COMMAND ${CMAKE_COMMAND} -E echo "The flatpak and flatpak-builder commands required to build a LAMMPS-GUI flatpak bundle were not found. Skipping.")
  endif()

  if(APPLE)
    file(STRINGS ${LAMMPS_DIR}/src/version.h line REGEX LAMMPS_VERSION)
    string(REGEX REPLACE "#define LAMMPS_VERSION \"([0-9]+) ([A-Za-z][A-Za-z][A-Za-z])[A-Za-z]* ([0-9]+)\""
                        "\\1\\2\\3" LAMMPS_RELEASE "${line}")

    # additional targets to populate the bundle tree and create the .dmg image file
    set(APP_CONTENTS ${CMAKE_BINARY_DIR}/lammps-gui_build-prefix/bin/lammps-gui.app/Contents)
    if(BUILD_TOOLS)
      file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/lammps-gui_build-prefix/bin/lammps-gui.app)
      add_custom_target(complete-bundle
        ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/bin
        COMMAND ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/Frameworks
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:lammps> ${APP_CONTENTS}/Frameworks/
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:lmp> ${APP_CONTENTS}/bin/
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/lmp ${APP_CONTENTS}/bin/
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/msi2lmp ${APP_CONTENTS}/bin/
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/binary2txt ${APP_CONTENTS}/bin/
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/stl_bin2txt ${APP_CONTENTS}/bin/
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/phana ${APP_CONTENTS}/bin/
        COMMAND ${CMAKE_COMMAND} -E create_symlink ../MacOS/lammps-gui ${APP_CONTENTS}/bin/lammps-gui
        COMMAND ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/Resources
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/cmake/packaging/README.macos ${APP_CONTENTS}/Resources/README.txt
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/cmake/packaging/lammps.icns ${APP_CONTENTS}/Resources
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/cmake/packaging/lammps-gui.icns ${APP_CONTENTS}/Resources
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/cmake/packaging/LAMMPS_DMG_Background.png ${APP_CONTENTS}/Resources
        COMMAND ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/share/lammps
        COMMAND ${CMAKE_COMMAND} -E make_directory ${APP_CONTENTS}/share/lammps/man/man1
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${LAMMPS_DIR}/potentials ${APP_CONTENTS}/share/lammps/potentials
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${LAMMPS_DIR}/bench ${APP_CONTENTS}/share/lammps/bench
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${LAMMPS_DIR}/tools/msi2lmp/frc_files ${APP_CONTENTS}/share/lammps/frc_files
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/doc/lammps.1 ${APP_CONTENTS}/share/lammps/man/man1/
        COMMAND ${CMAKE_COMMAND} -E create_symlink lammps.1 ${APP_CONTENTS}/share/lammps/man/man1/lmp.1
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${LAMMPS_DIR}/doc/msi2lmp.1 ${APP_CONTENTS}/share/lammps/man/man1
        DEPENDS lammps lmp tools lammps-gui_build
        COMMENT "Copying additional files into macOS app bundle tree"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      )
    else()
      message(FATAL_ERROR "Must use -D BUILD_TOOLS=yes for building app bundle")
    endif()
    if(BUILD_WHAM)
      add_custom_target(copy-wham
        ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/wham ${APP_CONTENTS}/bin/
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/wham-2d ${APP_CONTENTS}/bin/
        DEPENDS complete-bundle wham wham-2d
        COMMENT "Copying WHAM executables into macOS app bundle tree"
      )
      set(WHAM_TARGET copy-wham)
    endif()
    if(FFMPEG_EXECUTABLE)
      add_custom_target(copy-ffmpeg
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${FFMPEG_EXECUTABLE} ${APP_CONTENTS}/bin/
        COMMENT "Copying FFMpeg into macOS app bundle tree"
        DEPENDS complete-bundle
      )
      set(FFMPEG_TARGET copy-ffmpeg)
    endif()
    add_custom_target(dmg
      COMMAND ${LAMMPS_DIR}/cmake/packaging/build_macos_dmg.sh ${LAMMPS_RELEASE} ${CMAKE_BINARY_DIR}/lammps-gui_build-prefix/bin/lammps-gui.app
      DEPENDS complete-bundle ${WHAM_TARGET} ${FFMPEG_TARGET}
      COMMENT "Create Drag-n-Drop installer disk image from app bundle"
      BYPRODUCT LAMMPS-macOS-multiarch-GUI-${LAMMPS_RELEASE}.dmg
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    # settings or building on Windows with Visual Studio
  elseif(MSVC)
    file(STRINGS ${LAMMPS_DIR}/src/version.h line REGEX LAMMPS_VERSION)
    string(REGEX REPLACE "#define LAMMPS_VERSION \"([0-9]+) ([A-Za-z][A-Za-z][A-Za-z])[A-Za-z]* ([0-9]+)\""
                          "\\1\\2\\3" LAMMPS_RELEASE "${line}")
    #    install(FILES $<TARGET_RUNTIME_DLLS:lammps-gui> TYPE BIN)
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
    file(STRINGS ${LAMMPS_DIR}/src/version.h line REGEX LAMMPS_VERSION)
    string(REGEX REPLACE "#define LAMMPS_VERSION \"([0-9]+) ([A-Za-z][A-Za-z][A-Za-z])[A-Za-z]* ([0-9]+)\""
                          "\\1\\2\\3" LAMMPS_RELEASE "${line}")
    if(BUILD_SHARED_LIBS)
      install(FILES $<TARGET_RUNTIME_DLLS:lammps> TYPE BIN)
    endif()
    install(FILES $<TARGET_RUNTIME_DLLS:lmp> TYPE BIN)
    add_custom_target(zip
      COMMAND sh -vx ${LAMMPS_DIR}/cmake/packaging/build_windows_cross_zip.sh ${CMAKE_INSTALL_PREFIX} ${LAMMPS_RELEASE}
      DEPENDS lmp lammps-gui_build ${WHAM_EXE}
      COMMENT "Create zip file with windows binaries"
      BYPRODUCT LAMMPS-Win10-amd64-${LAMMPS_VERSION}.zip
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  elseif((CMAKE_SYSTEM_NAME STREQUAL "Linux") AND NOT LAMMPS_GUI_USE_PLUGIN)
    file(STRINGS ${LAMMPS_DIR}/src/version.h line REGEX LAMMPS_VERSION)
    string(REGEX REPLACE "#define LAMMPS_VERSION \"([0-9]+) ([A-Za-z][A-Za-z][A-Za-z])[A-Za-z]* ([0-9]+)\""
      "\\1\\2\\3" LAMMPS_RELEASE "${line}")
    set(LAMMPS_GUI_PACKAGING ${CMAKE_BINARY_DIR}/lammps-gui_build-prefix/src/lammps-gui_build/packaging/)
    set(LAMMPS_GUI_RESOURCES ${CMAKE_BINARY_DIR}/lammps-gui_build-prefix/src/lammps-gui_build/resources/)
    install(PROGRAMS ${CMAKE_BINARY_DIR}/lammps-gui_build-prefix/bin/lammps-gui DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(FILES ${LAMMPS_GUI_PACKAGING}/lammps-gui.desktop DESTINATION ${CMAKE_INSTALL_DATADIR}/applications/)
    install(FILES ${LAMMPS_GUI_PACKAGING}/lammps-gui.appdata.xml DESTINATION ${CMAKE_INSTALL_DATADIR}/appdata/)
    install(FILES ${LAMMPS_GUI_PACKAGING}/lammps-input.xml DESTINATION ${CMAKE_INSTALL_DATADIR}/mime/packages/)
    install(FILES ${LAMMPS_GUI_PACKAGING}/lammps-input.xml DESTINATION ${CMAKE_INSTALL_DATADIR}/mime/text/x-application-lammps.xml)
    install(DIRECTORY ${LAMMPS_GUI_RESOURCES}/icons/hicolor DESTINATION ${CMAKE_INSTALL_DATADIR}/icons/)
    install(CODE [[
      file(GET_RUNTIME_DEPENDENCIES
        LIBRARIES $<TARGET_FILE:lammps>
        EXECUTABLES $<TARGET_FILE:lmp> ${CMAKE_BINARY_DIR}/lammps-gui_build-prefix/bin/lammps-gui
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

    if(USE_INTERNAL_LINALG AND (NOT DOWNLOAD_POTENTIALS))
      add_custom_target(tgz
        COMMAND ${LAMMPS_DIR}/cmake/packaging/build_linux_tgz.sh ${LAMMPS_RELEASE}
        DEPENDS lmp tools lammps-gui_build ${WHAM_EXE}
        COMMENT "Create compressed tar file of LAMMPS-GUI with dependent libraries and wrapper"
        BYPRODUCT LAMMPS-Linux-x86_64-GUI-${LAMMPS_RELEASE}.tar.gz
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      )
    else()
      if(DOWNLOAD_POTENTIALS)
        add_custom_target(tgz
              COMMAND ${CMAKE_COMMAND} -E echo "Must use -D DOWLOAD_POTENTIALS=OFF for building Linux tgz package")
      else()
        add_custom_target(tgz
              COMMAND ${CMAKE_COMMAND} -E echo "Must use -D USE_INTERNAL_LINALG=ON for building Linux tgz package")
      endif()
    endif()
  endif()
endif()
