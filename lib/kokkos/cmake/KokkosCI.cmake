cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

message(STATUS "")

get_cmake_property(_cached_vars CACHE_VARIABLES)
set(KOKKOS_CMAKE_ARGS)
set(EXCLUDED_VARIABLES "CMAKE_COMMAND" "CMAKE_CPACK_COMMAND" "CMAKE_CTEST_COMMAND" "CMAKE_ROOT"
                       "CTEST_ARGS" "BUILD_NAME" "CMAKE_CXX_FLAGS" "CMAKE_BUILD_TYPE")
list(SORT _cached_vars)
foreach(_var ${_cached_vars})
    if(NOT "${_var}" IN_LIST EXCLUDED_VARIABLES)
        list(APPEND KOKKOS_CMAKE_ARGS ${_var})
        if("${_var}" STREQUAL "CMAKE_BUILD_TYPE")
            set(BUILD_TYPE "${CMAKE_BUILD_TYPE}")
        endif()
    endif()
endforeach()


#----------------------------------------------------------------------------------------#
#
#   Macros and variables
#
#----------------------------------------------------------------------------------------#

macro(CHECK_REQUIRED VAR)
    if(NOT DEFINED ${VAR})
        message(FATAL_ERROR "Error! Variable '${VAR}' must be defined")
    endif()
endmacro()

# require the build name variable
CHECK_REQUIRED(BUILD_NAME)

# uses all args
macro(SET_DEFAULT VAR)
    if(NOT DEFINED ${VAR})
        set(${VAR} ${ARGN})
    endif()
    # remove these ctest configuration variables from the defines
    # passed to the Kokkos configuration
    if("${VAR}" IN_LIST KOKKOS_CMAKE_ARGS)
        list(REMOVE_ITEM KOKKOS_CMAKE_ARGS "${VAR}")
    endif()
endmacro()

# uses first arg -- useful for selecting via priority from multiple
# potentially defined variables, e.g.:
#
#   set_default_arg1(BUILD_NAME ${TRAVIS_BUILD_NAME} ${BUILD_NAME})
#
macro(SET_DEFAULT_ARG1 VAR)
    if(NOT DEFINED ${VAR})
        foreach(_ARG ${ARGN})
            if(NOT "${_ARG}" STREQUAL "")
                set(${VAR} ${_ARG})
                break()
            endif()
        endforeach()
    endif()
    # remove these ctest configuration variables from the defines
    # passed to the Kokkos configuration
    if("${VAR}" IN_LIST KOKKOS_CMAKE_ARGS)
        list(REMOVE_ITEM KOKKOS_CMAKE_ARGS "${VAR}")
    endif()
endmacro()

# determine the default working directory
if(NOT "$ENV{WORKSPACE}" STREQUAL "")
    set(WORKING_DIR "$ENV{WORKSPACE}")
else()
    get_filename_component(WORKING_DIR ${CMAKE_CURRENT_LIST_DIR} DIRECTORY)
endif()

# determine the hostname
execute_process(COMMAND hostname
    OUTPUT_VARIABLE HOSTNAME
    OUTPUT_STRIP_TRAILING_WHITESPACE)

SET_DEFAULT(HOSTNAME "$ENV{HOSTNAME}")

# get the number of processors
include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)

# find git
find_package(Git QUIET)
if(NOT GIT_EXECUTABLE)
    unset(GIT_EXECUTABLE CACHE)
    unset(GIT_EXECUTABLE)
endif()

function(EXECUTE_GIT_COMMAND VAR)
    set(${VAR} "" PARENT_SCOPE)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${ARGN}
        OUTPUT_VARIABLE VAL
        RESULT_VARIABLE RET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        ERROR_QUIET)
    string(REPLACE ";" " " _CMD "${GIT_EXECUTABLE} ${ARGN}")
    set(LAST_GIT_COMMAND "${_CMD}" PARENT_SCOPE)
    if(RET EQUAL 0)
        set(${VAR} "${VAL}" PARENT_SCOPE)
    endif()
endfunction()

# just gets the git branch name if available
function(GET_GIT_BRANCH_NAME VAR)
    execute_git_command(GIT_BRANCH branch --show-current)
    set(_INVALID "%D" "HEAD")
    if(NOT GIT_BRANCH OR "${GIT_BRANCH}" IN_LIST _INVALID)
        execute_git_command(GIT_BRANCH show -s --format=%D)
        if(NOT GIT_BRANCH OR "${GIT_BRANCH}" IN_LIST _INVALID)
            execute_git_command(GIT_BRANCH --describe all)
        endif()
    endif()
    #
    if(GIT_BRANCH)
        string(REPLACE " " ";" _DESC "${GIT_BRANCH}")
        # just set it to last one via loop instead of wonky cmake index manip
        foreach(_ITR ${_DESC})
            set(GIT_BRANCH "${_ITR}")
        endforeach()
        set(${VAR} "${GIT_BRANCH}" PARENT_SCOPE)
        message(STATUS "GIT BRANCH via '${LAST_GIT_COMMAND}': ${GIT_BRANCH}")
    endif()
endfunction()

# just gets the git branch name if available
function(GET_GIT_AUTHOR_NAME VAR)
    execute_git_command(GIT_AUTHOR show -s --format=%an)
    if(GIT_AUTHOR)
        string(LENGTH "${GIT_AUTHOR}" STRLEN)
        # if the build name gets too long, this can cause submission errors
        if(STRLEN GREATER 24)
            # remove middle initial
            string(REGEX REPLACE " [A-Z]\. " " " GIT_AUTHOR "${GIT_AUTHOR}")
            # get first and sur name
            string(REGEX REPLACE "([A-Za-z]+) ([A-Za-z]+)" "\\1" F_NAME "${GIT_AUTHOR}")
            string(REGEX REPLACE "([A-Za-z]+) ([A-Za-z]+)" "\\2" S_NAME "${GIT_AUTHOR}")
            if(S_NAME)
                set(GIT_AUTHOR "${S_NAME}")
            elseif(F_NAME)
                set(GIT_AUTHOR "${F_NAME}")
            endif()
        endif()
        # remove any spaces, quotes, periods, etc.
        string(REGEX REPLACE "[ ',;_\.\"]+" "" GIT_AUTHOR "${GIT_AUTHOR}")
        set(${VAR} "${GIT_AUTHOR}" PARENT_SCOPE)
        message(STATUS "GIT AUTHOR via '${LAST_GIT_COMMAND}': ${GIT_AUTHOR}")
    endif()
endfunction()

# get the name of the branch
GET_GIT_BRANCH_NAME(GIT_BRANCH)
# get the name of the author
GET_GIT_AUTHOR_NAME(GIT_AUTHOR)
# author, prefer git method for consistency
SET_DEFAULT_ARG1(AUTHOR ${GIT_AUTHOR} $ENV{GIT_AUTHOR} $ENV{AUTHOR})
# SLUG == owner_name/repo_name
SET_DEFAULT_ARG1(SLUG $ENV{TRAVIS_PULL_REQUEST_SLUG} $ENV{TRAVIS_REPO_SLUG} $ENV{APPVEYOR_REPO_NAME} $ENV{PULL_REQUEST_SLUG} $ENV{REPO_SLUG})
# branch name
SET_DEFAULT_ARG1(BRANCH $ENV{TRAVIS_PULL_REQUEST_BRANCH} $ENV{TRAVIS_BRANCH} $ENV{APPVEYOR_PULL_REQUEST_HEAD_REPO_BRANCH} $ENV{APPVEYOR_REPO_BRANCH} $ENV{GIT_BRANCH} $ENV{BRANCH_NAME} $ENV{BRANCH} ${GIT_BRANCH})
# pull request number
SET_DEFAULT_ARG1(PULL_REQUEST_NUM $ENV{TRAVIS_PULL_REQUEST} $ENV{CHANGE_ID} $ENV{APPVEYOR_PULL_REQUEST_NUMBER} $ENV{PULL_REQUEST_NUM})
# get the event type, e.g. push, pull_request, api, cron, etc.
SET_DEFAULT_ARG1(EVENT_TYPE $ENV{TRAVIS_EVENT_TYPE} ${EVENT_TYPE})

if("${BRANCH}" STREQUAL "")
    message(STATUS "Checked: environment variables for Travis, Appveyor, Jenkins (git plugin), BRANCH_NAME, BRANCH and 'git branch --show-current'")
    message(FATAL_ERROR "Error! Git branch could not be determined. Please provide -DBRANCH=<name>")
endif()

#----------------------------------------------------------------------------------------#
#
#   Set default values if not provided on command-line
#
#----------------------------------------------------------------------------------------#

SET_DEFAULT(SOURCE_DIR      "${WORKING_DIR}")           # source directory
SET_DEFAULT(BINARY_DIR      "${WORKING_DIR}/build")     # build directory
SET_DEFAULT(BUILD_TYPE      "${CMAKE_BUILD_TYPE}")      # Release, Debug, etc.
SET_DEFAULT(MODEL           "Continuous")               # Continuous, Nightly, or Experimental
SET_DEFAULT(JOBS            1)                          # number of parallel ctests
SET_DEFAULT(CTEST_COMMAND   "${CMAKE_CTEST_COMMAND}")   # just in case
SET_DEFAULT(CTEST_ARGS      "-V --output-on-failure")   # extra arguments when ctest is called
SET_DEFAULT(GIT_EXECUTABLE  "git")                      # ctest_update
SET_DEFAULT(TARGET          "all")                      # build target
SET_DEFAULT_ARG1(SITE       "$ENV{SITE}"
                            "${HOSTNAME}")              # update site
SET_DEFAULT_ARG1(BUILD_JOBS "$ENV{BUILD_JOBS}"
                            "${NUM_PROCESSORS}")        # number of parallel compile jobs
#
#   The variable below correspond to ctest arguments, i.e. START,END,STRIDE are
#   '-I START,END,STRIDE'
#
SET_DEFAULT(START           "")
SET_DEFAULT(END             "")
SET_DEFAULT(STRIDE          "")
SET_DEFAULT(INCLUDE         "")
SET_DEFAULT(EXCLUDE         "")
SET_DEFAULT(INCLUDE_LABEL   "")
SET_DEFAULT(EXCLUDE_LABEL   "")
SET_DEFAULT(PARALLEL_LEVEL  "")
SET_DEFAULT(STOP_TIME       "")
SET_DEFAULT(LABELS          "")
SET_DEFAULT(NOTES           "")

# default static build tag for Nightly
set(BUILD_TAG "${BRANCH}")

if(NOT BUILD_TYPE)
    # default for kokkos if not specified
    set(BUILD_TYPE "RelWithDebInfo")
endif()

# generate dynamic name if continuous or experimental model
if(NOT "${MODEL}" STREQUAL "Nightly")
    if(EVENT_TYPE AND PULL_REQUEST_NUM)
        # e.g. pull_request/123
        if(AUTHOR)
            set(BUILD_TAG "${AUTHOR}/${EVENT_TYPE}/${PULL_REQUEST_NUM}")
        else()
            set(BUILD_TAG "${EVENT_TYPE}/${PULL_REQUEST_NUM}")
        endif()
    elseif(SLUG)
        # e.g. owner_name/repo_name
        set(BUILD_TAG "${SLUG}")
    elseif(AUTHOR)
        set(BUILD_TAG "${AUTHOR}/${BRANCH}")
    endif()
    if(EVENT_TYPE AND NOT PULL_REQUEST_NUM)
        set(BUILD_TAG "${BUILD_TAG}-${EVENT_TYPE}")
    endif()
endif()

# unnecessary
string(REPLACE "/remotes/" "/" BUILD_TAG "${BUILD_TAG}")
string(REPLACE "/origin/" "/" BUILD_TAG "${BUILD_TAG}")

message(STATUS "BUILD_TAG: ${BUILD_TAG}")

set(BUILD_NAME "[${BUILD_TAG}] [${BUILD_NAME}-${BUILD_TYPE}]")

# colons in build name create extra (empty) entries in CDash
string(REPLACE ":" "-" BUILD_NAME "${BUILD_NAME}")
# unnecessary info
string(REPLACE "/merge]" "]" BUILD_NAME "${BUILD_NAME}")
# consistency
string(REPLACE "/pr/" "/pull/" BUILD_NAME "${BUILD_NAME}")
string(REPLACE "pull_request/" "pull/" BUILD_NAME "${BUILD_NAME}")
# miscellaneous from missing fields
string(REPLACE "--" "-" BUILD_NAME "${BUILD_NAME}")
string(REPLACE "-]" "]" BUILD_NAME "${BUILD_NAME}")

# check binary directory
if(EXISTS ${BINARY_DIR})
    if(NOT IS_DIRECTORY "${BINARY_DIR}")
        message(FATAL_ERROR "Error! '${BINARY_DIR}' already exists and is not a directory!")
    endif()
    file(GLOB BINARY_DIR_FILES "${BINARY_DIR}/*")
    if(NOT "${BINARY_DIR_FILES}" STREQUAL "")
        message(FATAL_ERROR "Error! '${BINARY_DIR}' already exists and is not empty!")
    endif()
endif()

get_filename_component(SOURCE_REALDIR ${SOURCE_DIR} REALPATH)
get_filename_component(BINARY_REALDIR ${BINARY_DIR} REALPATH)

#----------------------------------------------------------------------------------------#
#
#   Generate the CTestConfig.cmake
#
#----------------------------------------------------------------------------------------#

set(CONFIG_ARGS)
foreach(_ARG ${KOKKOS_CMAKE_ARGS})
    if(NOT "${${_ARG}}" STREQUAL "")
        get_property(_ARG_TYPE CACHE ${_ARG} PROPERTY TYPE)
        if("${_ARG_TYPE}" STREQUAL "UNINITIALIZED")
            if("${${_ARG}}" STREQUAL "ON" OR "${${_ARG}}" STREQUAL "OFF")
                set(_ARG_TYPE "BOOL")
            elseif(EXISTS "${${_ARG}}" AND NOT IS_DIRECTORY "${${_ARG}}")
                set(_ARG_TYPE "FILEPATH")
            elseif(EXISTS "${${_ARG}}" AND IS_DIRECTORY "${${_ARG}}")
                set(_ARG_TYPE "PATH")
            elseif(NOT "${${_ARG}}" STREQUAL "")
                set(_ARG_TYPE "STRING")
            endif()
        endif()
        set(CONFIG_ARGS "${CONFIG_ARGS}set(${_ARG} \"${${_ARG}}\" CACHE ${_ARG_TYPE} \"\")\n")
    endif()
endforeach()

file(WRITE ${BINARY_REALDIR}/initial-cache.cmake
"
set(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS}\" CACHE STRING \"\")
${CONFIG_ARGS}
")

file(READ ${BINARY_REALDIR}/initial-cache.cmake _CACHE_INFO)
message(STATUS "Initial cache:\n${_CACHE_INFO}")

# initialize the cache
set(CONFIG_ARGS "-C ${BINARY_REALDIR}/initial-cache.cmake")


# generate the CTestConfig.cmake
configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/CTestConfig.cmake.in
    ${BINARY_REALDIR}/CTestConfig.cmake
    @ONLY)

# copy/generate the dashboard script
configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/KokkosCTest.cmake.in
    ${BINARY_REALDIR}/KokkosCTest.cmake
    @ONLY)

# custom CTest settings go in ${BINARY_DIR}/CTestCustom.cmake
execute_process(
    COMMAND             ${CMAKE_COMMAND} -E touch CTestCustom.cmake
    WORKING_DIRECTORY   ${BINARY_REALDIR}
    )

#----------------------------------------------------------------------------------------#
#
#   Execute CTest
#
#----------------------------------------------------------------------------------------#

message(STATUS "")
message(STATUS "BUILD_NAME: ${BUILD_NAME}")
message(STATUS "Executing '${CTEST_COMMAND} -S KokkosCTest.cmake ${CTEST_ARGS}'...")
message(STATUS "")

# e.g. -DCTEST_ARGS="--output-on-failure -VV" should really be -DCTEST_ARGS="--output-on-failure;-VV"
string(REPLACE " " ";" CTEST_ARGS "${CTEST_ARGS}")

execute_process(
    COMMAND             ${CTEST_COMMAND} -S KokkosCTest.cmake ${CTEST_ARGS}
    RESULT_VARIABLE     RET
    WORKING_DIRECTORY   ${BINARY_REALDIR}
    )

# ensure that any non-zero result variable gets propagated
if(NOT RET EQUAL 0)
    message(FATAL_ERROR "CTest return non-zero exit code: ${RET}")
endif()
