#These are tribits wrappers used by all projects in the Kokkos ecosystem

include(CMakeParseArguments)
include(CTest)

function(ASSERT_DEFINED VARS)
  foreach(VAR ${VARS})
    if(NOT DEFINED ${VAR})
      message(SEND_ERROR "Error, the variable ${VAR} is not defined!")
    endif()
  endforeach()
endfunction()

macro(APPEND_GLOB VAR)
  file(GLOB LOCAL_TMP_VAR ${ARGN})
  list(APPEND ${VAR} ${LOCAL_TMP_VAR})
endmacro()

macro(GLOBAL_SET VARNAME)
  set(${VARNAME} ${ARGN} CACHE INTERNAL "" FORCE)
endmacro()

macro(PREPEND_GLOBAL_SET VARNAME)
  assert_defined(${VARNAME})
  global_set(${VARNAME} ${ARGN} ${${VARNAME}})
endmacro()

macro(ADD_INTERFACE_LIBRARY LIB_NAME)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp "")
  add_library(${LIB_NAME} STATIC ${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp)
  set_target_properties(${LIB_NAME} PROPERTIES INTERFACE TRUE)
endmacro()

function(KOKKOS_ADD_TEST)
  cmake_parse_arguments(
    TEST "WILL_FAIL;SKIP_TRIBITS" "FAIL_REGULAR_EXPRESSION;PASS_REGULAR_EXPRESSION;EXE;NAME;TOOL" "CATEGORIES;ARGS"
    ${ARGN}
  )
  # To match Tribits, we should always be receiving
  # the root names of exes/libs
  if(TEST_EXE)
    set(EXE_ROOT ${TEST_EXE})
  else()
    set(EXE_ROOT ${TEST_NAME})
  endif()
  # Prepend package name to the test name
  # These should be the full target name
  set(TEST_NAME ${PACKAGE_NAME}_${TEST_NAME})

  # For compatibility with Trilinos testing, we support:
  #  * `-D <fullTestName>_DISABLE=ON`
  #  * `-D <fullTestName>_EXTRA_ARGS="<arg0>;<arg1>;<arg2>;..."`
  #  * `-D <fullTestName>_SET_RUN_SERIAL=ON`
  if(${TEST_NAME}_DISABLE)
    return()
  endif()

  set(EXE ${PACKAGE_NAME}_${EXE_ROOT})
  if(WIN32)
    add_test(NAME ${TEST_NAME} WORKING_DIRECTORY ${LIBRARY_OUTPUT_PATH} COMMAND ${EXE}${CMAKE_EXECUTABLE_SUFFIX}
                                                                                ${TEST_ARGS} ${${TEST_NAME}_EXTRA_ARGS}
    )
  else()
    add_test(NAME ${TEST_NAME} COMMAND ${EXE} ${TEST_ARGS} ${${TEST_NAME}_EXTRA_ARGS})
  endif()
  # Trilinos testing benefits from labeling the tests as "Kokkos" tests
  set_tests_properties(${TEST_NAME} PROPERTIES LABELS Kokkos)
  if(${TEST_NAME}_SET_RUN_SERIAL)
    set_tests_properties(${TEST_NAME} PROPERTIES RUN_SERIAL ON)
  endif()
  # TriBITS doesn't actually currently support `-D <fullTestName>_ENVIRONMENT`
  # but we decided to add it anyway
  if(${TEST_NAME}_ENVIRONMENT)
    set_tests_properties(${TEST_NAME} PROPERTIES ENVIRONMENT "${${TEST_NAME}_ENVIRONMENT}")
  endif()
  if(TEST_WILL_FAIL)
    set_tests_properties(${TEST_NAME} PROPERTIES WILL_FAIL ${TEST_WILL_FAIL})
  endif()
  if(TEST_FAIL_REGULAR_EXPRESSION)
    set_tests_properties(${TEST_NAME} PROPERTIES FAIL_REGULAR_EXPRESSION ${TEST_FAIL_REGULAR_EXPRESSION})
  endif()
  if(TEST_PASS_REGULAR_EXPRESSION)
    set_tests_properties(${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION ${TEST_PASS_REGULAR_EXPRESSION})
  endif()
  if(TEST_TOOL)
    add_dependencies(${EXE} ${TEST_TOOL}) #make sure the exe has to build the tool
    set_property(TEST ${TEST_NAME} APPEND_STRING PROPERTY ENVIRONMENT "KOKKOS_TOOLS_LIBS=$<TARGET_FILE:${TEST_TOOL}>")
  endif()
  verify_empty(KOKKOS_ADD_TEST ${TEST_UNPARSED_ARGUMENTS})
endfunction()

macro(KOKKOS_CREATE_IMPORTED_TPL_LIBRARY TPL_NAME)
  add_interface_library(TPL_LIB_${TPL_NAME})
  target_link_libraries(TPL_LIB_${TPL_NAME} LINK_PUBLIC ${TPL_${TPL_NAME}_LIBRARIES})
  target_include_directories(TPL_LIB_${TPL_NAME} INTERFACE ${TPL_${TPL_NAME}_INCLUDE_DIRS})
endmacro()

function(KOKKOS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES TPL_NAME)
  cmake_parse_arguments(PARSE "" "" "REQUIRED_HEADERS;REQUIRED_LIBS_NAMES" ${ARGN})

  set(_${TPL_NAME}_ENABLE_SUCCESS TRUE)
  if(PARSE_REQUIRED_LIBS_NAMES)
    find_library(TPL_${TPL_NAME}_LIBRARIES NAMES ${PARSE_REQUIRED_LIBS_NAMES})
    if(NOT TPL_${TPL_NAME}_LIBRARIES)
      set(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
    endif()
  endif()
  if(PARSE_REQUIRED_HEADERS)
    find_path(TPL_${TPL_NAME}_INCLUDE_DIRS NAMES ${PARSE_REQUIRED_HEADERS})
    if(NOT TPL_${TPL_NAME}_INCLUDE_DIRS)
      set(_${TPL_NAME}_ENABLE_SUCCESS FALSE)
    endif()
  endif()
  if(_${TPL_NAME}_ENABLE_SUCCESS)
    kokkos_create_imported_tpl_library(${TPL_NAME})
  endif()
  verify_empty(KOKKOS_CREATE_IMPORTED_TPL_LIBRARY ${PARSE_UNPARSED_ARGUMENTS})
endfunction()

function(KOKKOS_LIB_TYPE LIB RET)
  get_target_property(PROP ${LIB} TYPE)
  if(${PROP} STREQUAL "INTERFACE_LIBRARY")
    set(${RET} "INTERFACE" PARENT_SCOPE)
  else()
    set(${RET} "PUBLIC" PARENT_SCOPE)
  endif()
endfunction()

function(KOKKOS_TARGET_INCLUDE_DIRECTORIES TARGET)
  if(TARGET ${TARGET})
    #the target actually exists - this means we are doing separate libs
    #or this a test library
    kokkos_lib_type(${TARGET} INCTYPE)
    target_include_directories(${TARGET} ${INCTYPE} ${ARGN})
  else()
    get_property(LIBS GLOBAL PROPERTY KOKKOS_LIBRARIES_NAMES)
    if(${TARGET} IN_LIST LIBS)
      set_property(GLOBAL APPEND PROPERTY KOKKOS_LIBRARY_INCLUDES ${ARGN})
    else()
      message(FATAL_ERROR "Trying to set include directories on unknown target ${TARGET}")
    endif()
  endif()
endfunction()

function(KOKKOS_LINK_INTERNAL_LIBRARY TARGET DEPLIB)
  set(options INTERFACE)
  set(oneValueArgs)
  set(multiValueArgs)
  cmake_parse_arguments(PARSE "INTERFACE" "" "" ${ARGN})
  set(LINK_TYPE)
  if(PARSE_INTERFACE)
    set(LINK_TYPE INTERFACE)
  else()
    set(LINK_TYPE PUBLIC)
  endif()
  target_link_libraries(${TARGET} ${LINK_TYPE} ${DEPLIB})
  verify_empty(KOKKOS_LINK_INTERNAL_LIBRARY ${PARSE_UNPARSED_ARGUMENTS})
endfunction()

function(KOKKOS_ADD_TEST_LIBRARY NAME)
  set(oneValueArgs)
  set(multiValueArgs HEADERS SOURCES)

  cmake_parse_arguments(PARSE "STATIC;SHARED" "" "HEADERS;SOURCES;DEPLIBS" ${ARGN})

  set(LIB_TYPE)
  if(PARSE_STATIC)
    set(LIB_TYPE STATIC)
  elseif(PARSE_SHARED)
    set(LIB_TYPE SHARED)
  endif()

  if(PARSE_HEADERS)
    list(REMOVE_DUPLICATES PARSE_HEADERS)
  endif()
  if(PARSE_SOURCES)
    list(REMOVE_DUPLICATES PARSE_SOURCES)
  endif()
  add_library(${NAME} ${LIB_TYPE} ${PARSE_SOURCES})
  if(PARSE_DEPLIBS)
    target_link_libraries(${NAME} PRIVATE ${PARSE_DEPLIBS})
  endif()
endfunction()

function(KOKKOS_INCLUDE_DIRECTORIES)
  cmake_parse_arguments(INC "REQUIRED_DURING_INSTALLATION_TESTING" "" "" ${ARGN})
  include_directories(${INC_UNPARSED_ARGUMENTS})
endfunction()

macro(PRINTALL match)
  get_cmake_property(_variableNames VARIABLES)
  list(SORT _variableNames)
  foreach(_variableName ${_variableNames})
    if("${_variableName}" MATCHES "${match}")
      message(STATUS "${_variableName}=${${_variableName}}")
    endif()
  endforeach()
endmacro()

macro(SET_GLOBAL_REPLACE SUBSTR VARNAME)
  string(REPLACE ${SUBSTR} ${${VARNAME}} TEMP)
  global_set(${VARNAME} ${TEMP})
endmacro()

function(GLOBAL_APPEND VARNAME)
  #We make this a function since we are setting variables
  #and want to use scope to avoid overwriting local variables
  set(TEMP ${${VARNAME}})
  list(APPEND TEMP ${ARGN})
  global_set(${VARNAME} ${TEMP})
endfunction()
