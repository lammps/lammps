include(FetchContent)

# PACE library support for ML-PACE package
find_package(pace QUIET)

if(pace_FOUND)
    target_link_libraries(lammps PRIVATE pace::pace)
else()
    # Silence FetchContent_Populate deprecation warning for CMake 4.0+
    if(POLICY CMP0169)
      cmake_policy(SET CMP0169 OLD)
    endif()
    # Silence warnings about timestamps of downloaded files
    if(POLICY CMP0135)
      cmake_policy(SET CMP0135 OLD)
    endif()

    # Downloading the specific branch: alphataubio-fitsnap
    FetchContent_Declare(
      pacelib
      GIT_REPOSITORY https://github.com/alphataubio/lammps-pyace.git
      GIT_TAG        alphataubio-fitsnap
      GIT_SHALLOW    OFF
    )

    if(LOCAL_ML-PACE)
      set(lib-pace "${LOCAL_ML-PACE}")
    else()
      # This clones the repo and checks out the branch
      FetchContent_Populate(pacelib)
      set(lib-pace "${pacelib_SOURCE_DIR}")
    endif()

    # Ensure yaml-cpp namespaced target exists
    find_package(yaml-cpp QUIET)
    if(TARGET yaml-cpp AND NOT TARGET yaml-cpp::yaml-cpp)
      add_library(yaml-cpp::yaml-cpp ALIAS yaml-cpp)
    endif()

    # Fixup yaml-cpp/emitterutils.cpp for GCC 15+ / AppleClang 17 compatibility
    if(EXISTS "${lib-pace}/yaml-cpp/src/emitterutils.cpp")
      file(READ ${lib-pace}/yaml-cpp/src/emitterutils.cpp yaml_emitterutils)
      string(REPLACE "#include <sstream>" "#include <sstream>\n#include <cinttypes>" yaml_tmp_emitterutils "${yaml_emitterutils}")
      string(REPLACE "#include <cinttypes>\n#include <cinttypes>" "#include <cinttypes>" yaml_emitterutils "${yaml_tmp_emitterutils}")
      file(WRITE ${lib-pace}/yaml-cpp/src/emitterutils.cpp "${yaml_emitterutils}")
    endif()

    add_subdirectory(${lib-pace} build-pace EXCLUDE_FROM_ALL)
    set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})

    if(CMAKE_PROJECT_NAME STREQUAL "lammps")
      target_link_libraries(lammps PRIVATE pace)
    endif()
endif()
