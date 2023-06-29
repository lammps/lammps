# https://jonathanhamberg.com/post/cmake-embedding-git-hash/

find_package(Git QUIET)

SET(CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_DIR})
SET(pre_configure_dir ${CMAKE_CURRENT_LIST_DIR})
SET(post_configure_dir ${CMAKE_BINARY_DIR}/generated)

SET(pre_configure_file ${pre_configure_dir}/Kokkos_Version_Info.cpp.in)
SET(post_configure_file ${post_configure_dir}/Kokkos_Version_Info.cpp)

FUNCTION(check_git_write git_hash git_clean_status)
  FILE(
    WRITE
    ${CMAKE_BINARY_DIR}/git-state.txt
    "${git_hash}-${git_clean_status}")
ENDFUNCTION()

FUNCTION(check_git_read git_hash)
  IF(EXISTS ${CMAKE_BINARY_DIR}/git-state.txt)
    FILE(STRINGS ${CMAKE_BINARY_DIR}/git-state.txt CONTENT)
    LIST(GET CONTENT 0 var)

    message(DEBUG "Cached Git hash: ${var}")
    SET(${git_hash} ${var} PARENT_SCOPE)
  else()
    SET(${git_hash} "INVALID" PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

FUNCTION(check_git_version)
  IF(NOT EXISTS ${post_configure_dir}/Kokkos_Version_Info.hpp)
    FILE(
      COPY ${pre_configure_dir}/Kokkos_Version_Info.hpp
      DESTINATION ${post_configure_dir})
  ENDIF()

  IF(NOT Git_FOUND OR NOT EXISTS ${KOKKOS_SOURCE_DIR}/.git)
    configure_file(${pre_configure_file} ${post_configure_file} @ONLY)
    return()
  ENDIF()

  # Get the current working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${KOKKOS_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Get the latest commit description
  execute_process(
    COMMAND ${GIT_EXECUTABLE} show -s --format=%s
    WORKING_DIRECTORY ${KOKKOS_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_DESCRIPTION
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Get the latest commit date
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%cI
    WORKING_DIRECTORY ${KOKKOS_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_DATE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Check if repo is dirty / clean
  execute_process(
    COMMAND ${GIT_EXECUTABLE} diff-index --quiet HEAD --
    WORKING_DIRECTORY ${KOKKOS_SOURCE_DIR}
    RESULT_VARIABLE IS_DIRTY
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  IF(IS_DIRTY EQUAL 0)
    SET(GIT_CLEAN_STATUS "CLEAN")
  else()
    SET(GIT_CLEAN_STATUS "DIRTY")
  ENDIF()

  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
    WORKING_DIRECTORY ${KOKKOS_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  check_git_read(GIT_HASH_CACHE)

  IF(NOT EXISTS ${post_configure_dir})
    file(MAKE_DIRECTORY ${post_configure_dir})
  ENDIF()

  # Only update the git_version.cpp if the hash has changed. This will
  # prevent us from rebuilding the project more than we need to.
  IF(NOT "${GIT_COMMIT_HASH}-${GIT_CLEAN_STATUS}" STREQUAL ${GIT_HASH_CACHE}
    OR NOT EXISTS ${post_configure_file})
    # Set the GIT_HASH_CACHE variable so the next build won't have
    # to regenerate the source file.
    check_git_write(${GIT_COMMIT_HASH} ${GIT_CLEAN_STATUS})

    configure_file(${pre_configure_file} ${post_configure_file} @ONLY)
    message(STATUS "Configured git information in ${post_configure_file}")
  ENDIF()
ENDFUNCTION()

FUNCTION(check_git_setup)
  add_custom_target(
    AlwaysCheckGit COMMAND ${CMAKE_COMMAND}
    -DRUN_CHECK_GIT_VERSION=1
    -DKOKKOS_SOURCE_DIR=${Kokkos_SOURCE_DIR}
    -P ${CURRENT_LIST_DIR}/build_env_info.cmake
    BYPRODUCTS ${post_configure_file})

  add_library(impl_git_version ${CMAKE_BINARY_DIR}/generated/Kokkos_Version_Info.cpp)
  target_include_directories(impl_git_version PUBLIC ${CMAKE_BINARY_DIR}/generated)
  target_compile_features(impl_git_version PRIVATE cxx_raw_string_literals)
  add_dependencies(impl_git_version AlwaysCheckGit)

  check_git_version()
ENDFUNCTION()

# This is used to run this function from an external cmake process.
IF(RUN_CHECK_GIT_VERSION)
  check_git_version()
ENDIF()
