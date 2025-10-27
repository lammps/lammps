################################### FUNCTIONS ##################################
# List of functions
#   kokkos_option

# Validate options are given with correct case and define an internal
# upper-case version for use within

set(Kokkos_OPTIONS_NOT_TO_EXPORT Kokkos_ENABLE_BENCHMARKS Kokkos_ENABLE_EXAMPLES Kokkos_ENABLE_TESTS
                                 Kokkos_ENABLE_HEADER_SELF_CONTAINMENT_TESTS Kokkos_ENABLE_COMPILER_WARNINGS
)

#
#
# @FUNCTION: kokkos_deprecated_list
#
# Function that checks if a deprecated list option like Kokkos_ARCH was given.
# This prints an error and prevents configure from completing.
# It attempts to print a helpful message about updating the options for the new CMake.
# Kokkos_${SUFFIX} is the name of the option (like Kokkos_ARCH) being checked.
# Kokkos_${PREFIX}_X is the name of new option to be defined from a list X,Y,Z,...
function(kokkos_deprecated_list SUFFIX PREFIX)
  set(CAMEL_NAME Kokkos_${SUFFIX})
  string(TOUPPER ${CAMEL_NAME} UC_NAME)

  #I don't love doing it this way but better to be safe
  foreach(opt ${KOKKOS_GIVEN_VARIABLES})
    string(TOUPPER ${opt} OPT_UC)
    if("${OPT_UC}" STREQUAL "${UC_NAME}")
      string(REPLACE "," ";" optlist "${${opt}}")
      set(ERROR_MSG
          "Given deprecated option list ${opt}. This must now be given as separate -D options, which assuming you spelled options correctly would be:"
      )
      foreach(entry ${optlist})
        string(TOUPPER ${entry} ENTRY_UC)
        string(APPEND ERROR_MSG "\n  -DKokkos_${PREFIX}_${ENTRY_UC}=ON")
      endforeach()
      string(
        APPEND
        ERROR_MSG
        "\nRemove CMakeCache.txt and re-run. For a list of valid options, refer to the configuration guide in the documention or even look at CMakeCache.txt (before deleting it)."
      )
      message(SEND_ERROR ${ERROR_MSG})
    endif()
  endforeach()
endfunction()

function(kokkos_option CAMEL_SUFFIX DEFAULT TYPE DOCSTRING)
  set(CAMEL_NAME Kokkos_${CAMEL_SUFFIX})
  string(TOUPPER ${CAMEL_NAME} UC_NAME)

  list(APPEND KOKKOS_OPTION_KEYS ${CAMEL_SUFFIX})
  set(KOKKOS_OPTION_KEYS ${KOKKOS_OPTION_KEYS} PARENT_SCOPE)
  list(APPEND KOKKOS_OPTION_VALUES "${DOCSTRING}")
  set(KOKKOS_OPTION_VALUES ${KOKKOS_OPTION_VALUES} PARENT_SCOPE)
  list(APPEND KOKKOS_OPTION_TYPES ${TYPE})
  set(KOKKOS_OPTION_TYPES ${KOKKOS_OPTION_TYPES} PARENT_SCOPE)

  # Make sure this appears in the cache with the appropriate DOCSTRING
  set(${CAMEL_NAME} ${DEFAULT} CACHE ${TYPE} ${DOCSTRING})

  #I don't love doing it this way because it's N^2 in number options, but c'est la vie
  foreach(opt ${KOKKOS_GIVEN_VARIABLES})
    string(TOUPPER ${opt} OPT_UC)
    if("${OPT_UC}" STREQUAL "${UC_NAME}")
      if(NOT "${opt}" STREQUAL "${CAMEL_NAME}")
        message(
          FATAL_ERROR
            "Matching option found for ${CAMEL_NAME} with the wrong case ${opt}. Please delete your CMakeCache.txt and change option to -D${CAMEL_NAME}=${${opt}}. This is now enforced to avoid hard-to-debug CMake cache inconsistencies."
        )
      endif()
    endif()
  endforeach()

  #okay, great, we passed the validation test - use the default
  if(DEFINED ${CAMEL_NAME})
    set(${UC_NAME} ${${CAMEL_NAME}} PARENT_SCOPE)
  else()
    set(${UC_NAME} ${DEFAULT} PARENT_SCOPE)
  endif()
endfunction()

include(CMakeDependentOption)
function(kokkos_dependent_option CAMEL_SUFFIX DOCSTRING DEFAULT DEPENDENCY FORCE)
  set(CAMEL_NAME Kokkos_${CAMEL_SUFFIX})
  string(TOUPPER ${CAMEL_NAME} UC_NAME)

  list(APPEND KOKKOS_OPTION_KEYS ${CAMEL_SUFFIX})
  set(KOKKOS_OPTION_KEYS ${KOKKOS_OPTION_KEYS} PARENT_SCOPE)
  list(APPEND KOKKOS_OPTION_VALUES "${DOCSTRING}")
  set(KOKKOS_OPTION_VALUES ${KOKKOS_OPTION_VALUES} PARENT_SCOPE)
  list(APPEND KOKKOS_OPTION_TYPES BOOL)
  set(KOKKOS_OPTION_TYPES ${KOKKOS_OPTION_TYPES} PARENT_SCOPE)

  cmake_dependent_option(${CAMEL_NAME} ${DOCSTRING} ${DEFAULT} "${DEPENDENCY}" ${FORCE})

  #I don't love doing it this way because it's N^2 in number options, but c'est la vie
  foreach(opt ${KOKKOS_GIVEN_VARIABLES})
    string(TOUPPER ${opt} OPT_UC)
    if("${OPT_UC}" STREQUAL "${UC_NAME}")
      if(NOT "${opt}" STREQUAL "${CAMEL_NAME}")
        message(
          FATAL_ERROR
            "Matching option found for ${CAMEL_NAME} with the wrong case ${opt}. Please delete your CMakeCache.txt and change option to -D${CAMEL_NAME}=${${opt}}. This is now enforced to avoid hard-to-debug CMake cache inconsistencies."
        )
      endif()
    endif()
  endforeach()

  #okay, great, we passed the validation test - use the default
  if(DEFINED ${CAMEL_NAME})
    set(${UC_NAME} ${${CAMEL_NAME}} PARENT_SCOPE)
  else()
    set(${UC_NAME} ${DEFAULT} PARENT_SCOPE)
  endif()
endfunction()

function(kokkos_set_option CAMEL_SUFFIX VALUE)
  list(FIND KOKKOS_OPTION_KEYS ${CAMEL_SUFFIX} OPTION_INDEX)
  if(OPTION_INDEX EQUAL -1)
    message(FATAL_ERROR "Couldn't set value for Kokkos_${CAMEL_SUFFIX}")
  endif()
  set(CAMEL_NAME Kokkos_${CAMEL_SUFFIX})
  string(TOUPPER ${CAMEL_NAME} UC_NAME)

  list(GET KOKKOS_OPTION_VALUES ${OPTION_INDEX} DOCSTRING)
  list(GET KOKKOS_OPTION_TYPES ${OPTION_INDEX} TYPE)
  set(${CAMEL_NAME} ${VALUE} CACHE ${TYPE} ${DOCSTRING} FORCE)
  message(STATUS "Setting ${CAMEL_NAME}=${VALUE}")
  set(${UC_NAME} ${VALUE} PARENT_SCOPE)
endfunction()

function(kokkos_append_config_line LINE)
  global_append(KOKKOS_TPL_EXPORTS "${LINE}")
endfunction()

macro(kokkos_export_cmake_tpl NAME)
  cmake_parse_arguments(KOKKOS_EXTRA_ARG "REQUIRED" "" "COMPONENTS" ${ARGN})

  #CMake TPLs are located with a call to find_package
  #find_package locates XConfig.cmake files through
  #X_DIR or X_ROOT variables set prior to calling find_package

  #If Kokkos was configured to find the TPL through a _DIR variable
  #make sure thar DIR variable is available to downstream packages
  if(DEFINED ${NAME}_DIR)
    #The downstream project may override the TPL location that Kokkos used
    #Check if the downstream project chose its own TPL location
    #If not, make the Kokkos found location available
    kokkos_append_config_line("IF(NOT DEFINED ${NAME}_DIR)")
    kokkos_append_config_line("  SET(${NAME}_DIR  ${${NAME}_DIR})")
    kokkos_append_config_line("ENDIF()")
  endif()

  if(DEFINED ${NAME}_ROOT)
    #The downstream project may override the TPL location that Kokkos used
    #Check if the downstream project chose its own TPL location
    #If not, make the Kokkos found location available
    kokkos_append_config_line("IF(NOT DEFINED ${NAME}_ROOT)")
    kokkos_append_config_line("  SET(${NAME}_ROOT  ${${NAME}_ROOT})")
    kokkos_append_config_line("ENDIF()")
  endif()
  set(KOKKOS_CONFIG_STRING "FIND_DEPENDENCY(${NAME}")

  if(KOKKOS_EXTRA_ARG_REQUIRED)
    string(APPEND KOKKOS_CONFIG_STRING " REQUIRED")
  endif()
  if(KOKKOS_EXTRA_ARG_COMPONENTS)
    string(APPEND KOKKOS_CONFIG_STRING " COMPONENTS ${KOKKOS_EXTRA_ARG_COMPONENTS}")
  endif()
  string(APPEND KOKKOS_CONFIG_STRING ")")
  kokkos_append_config_line(${KOKKOS_CONFIG_STRING})
endmacro()

macro(kokkos_export_imported_tpl NAME)
  get_target_property(LIB_IMPORTED ${NAME} IMPORTED)
  if(NOT LIB_IMPORTED)
    # This is not an imported target
    # This an interface library that we created
    install(
      TARGETS ${NAME}
      EXPORT KokkosTargets
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  else()
    #make sure this also gets "exported" in the config file
    kokkos_append_config_line("IF(NOT TARGET ${NAME})")

    get_target_property(LIB_TYPE ${NAME} TYPE)
    if(${LIB_TYPE} STREQUAL "INTERFACE_LIBRARY")
      kokkos_append_config_line("ADD_LIBRARY(${NAME} INTERFACE IMPORTED)")
      kokkos_append_config_line("SET_TARGET_PROPERTIES(${NAME} PROPERTIES")
    else()
      kokkos_append_config_line("ADD_LIBRARY(${NAME} UNKNOWN IMPORTED)")
      kokkos_append_config_line("SET_TARGET_PROPERTIES(${NAME} PROPERTIES")
      get_target_property(TPL_LIBRARY ${NAME} IMPORTED_LOCATION)
      if(TPL_LIBRARY)
        kokkos_append_config_line("IMPORTED_LOCATION \"${TPL_LIBRARY}\"")
      endif()
    endif()

    get_target_property(TPL_INCLUDES ${NAME} INTERFACE_INCLUDE_DIRECTORIES)
    if(TPL_INCLUDES)
      kokkos_append_config_line("INTERFACE_INCLUDE_DIRECTORIES \"${TPL_INCLUDES}\"")
    endif()

    get_target_property(TPL_COMPILE_OPTIONS ${NAME} INTERFACE_COMPILE_OPTIONS)
    if(TPL_COMPILE_OPTIONS)
      kokkos_append_config_line("INTERFACE_COMPILE_OPTIONS ${TPL_COMPILE_OPTIONS}")
    endif()

    set(TPL_LINK_OPTIONS)
    get_target_property(TPL_LINK_OPTIONS ${NAME} INTERFACE_LINK_OPTIONS)
    if(TPL_LINK_OPTIONS)
      kokkos_append_config_line("INTERFACE_LINK_OPTIONS ${TPL_LINK_OPTIONS}")
    endif()

    get_target_property(TPL_LINK_LIBRARIES ${NAME} INTERFACE_LINK_LIBRARIES)
    if(TPL_LINK_LIBRARIES)
      kokkos_append_config_line("INTERFACE_LINK_LIBRARIES \"${TPL_LINK_LIBRARIES}\"")
    endif()
    kokkos_append_config_line(")")
    kokkos_append_config_line("ENDIF()")
  endif()
endmacro()

#
# @MACRO: KOKKOS_IMPORT_TPL()
#
# Function that checks if a third-party library (TPL) has been enabled and calls `find_package`
# to create an imported target encapsulating all the flags and libraries
# needed to use the TPL
#
# Usage::
#
#   KOKKOS_IMPORT_TPL(
#     <NAME>
#     NO_EXPORT
#     INTERFACE
#
#   ``NO_EXPORT``
#
#     If specified, this TPL will not be added to KokkosConfig.cmake as an export
#
#   ``INTERFACE``
#
#     If specified, this TPL will build an INTERFACE library rather than an
#     IMPORTED target
macro(kokkos_import_tpl NAME)
  cmake_parse_arguments(TPL "NO_EXPORT;INTERFACE" "" "" ${ARGN})
  if(TPL_INTERFACE)
    set(TPL_IMPORTED_NAME ${NAME})
  else()
    set(TPL_IMPORTED_NAME Kokkos::${NAME})
  endif()

  if(KOKKOS_ENABLE_${NAME})
    #Tack on a TPL here to make sure we avoid using anyone else's find
    find_package(TPL${NAME} REQUIRED MODULE)
    if(NOT TARGET ${TPL_IMPORTED_NAME})
      message(FATAL_ERROR "Find module succeeded for ${NAME}, but did not produce valid target ${TPL_IMPORTED_NAME}")
    endif()
    if(NOT TPL_NO_EXPORT)
      get_target_property(TPL_ORIGINAL_NAME ${TPL_IMPORTED_NAME} ALIASED_TARGET)
      if(NOT TPL_ORIGINAL_NAME)
        set(TPL_ORIGINAL_NAME ${TPL_IMPORTED_NAME})
      endif()
      kokkos_export_imported_tpl(${TPL_ORIGINAL_NAME})
    endif()
    list(APPEND KOKKOS_ENABLED_TPLS ${NAME})
  endif()
endmacro(kokkos_import_tpl)

macro(kokkos_import_cmake_tpl MODULE_NAME)
  kokkos_import_tpl(${MODULE_NAME} ${ARGN} NO_EXPORT)
  cmake_parse_arguments(TPL "NO_EXPORT" "OPTION_NAME" "" ${ARGN})

  if(NOT TPL_OPTION_NAME)
    set(TPL_OPTION_NAME ${MODULE_NAME})
  endif()

  if(NOT TPL_NO_EXPORT)
    kokkos_export_cmake_tpl(${MODULE_NAME})
  endif()
endmacro()

#
# @MACRO: KOKKOS_CREATE_IMPORTED_TPL()
#
# Function that creates an imported target encapsulating all the flags
# and libraries needed to use the TPL
#
# Usage::
#
#   KOKKOS_CREATE_IMPORTED_TPL(
#     <NAME>
#     INTERFACE
#     LIBRARY <path_to_librarY>
#     LINK_LIBRARIES <lib1> <lib2> ...
#     COMPILE_OPTIONS <opt1> <opt2> ...
#     LINK_OPTIONS <opt1> <opt2> ...
#
#   ``INTERFACE``
#
#     If specified, this TPL will build an INTERFACE library rather than an
#     IMPORTED target
#
#   ``LIBRARY <path_to_library>``
#
#     If specified, this gives the IMPORTED_LOCATION of the library.
#
#   ``LINK_LIBRARIES <lib1> <lib2> ...``
#
#     If specified, this gives a list of dependent libraries that also
#     need to be linked against. Each entry can be a library path or
#     the name of a valid CMake target.
#
#   ``INCLUDES <path1> <path2> ...``
#
#     If specified, this gives a list of directories that must be added
#     to the include path for using this library.
#
#   ``COMPILE_OPTIONS <opt1> <opt2> ...``
#
#     If specified, this gives a list of compiler flags that must be used
#     for using this library.
#
#   ``LINK_OPTIONS <opt1> <opt2> ...``
#
#     If specified, this gives a list of linker flags that must be used
#     for using this library.
macro(kokkos_create_imported_tpl NAME)
  cmake_parse_arguments(
    TPL "INTERFACE" "LIBRARY" "LINK_LIBRARIES;INCLUDES;COMPILE_DEFINITIONS;COMPILE_OPTIONS;LINK_OPTIONS" ${ARGN}
  )

  if(TPL_INTERFACE)
    add_library(${NAME} INTERFACE)
    #Give this an importy-looking name
    add_library(Kokkos::${NAME} ALIAS ${NAME})
    if(TPL_LIBRARY)
      message(SEND_ERROR "TPL Interface library ${NAME} should not have an IMPORTED_LOCATION")
    endif()
    #Things have to go in quoted in case we have multiple list entries
    if(TPL_LINK_LIBRARIES)
      target_link_libraries(${NAME} INTERFACE ${TPL_LINK_LIBRARIES})
    endif()
    if(TPL_INCLUDES)
      target_include_directories(${NAME} INTERFACE ${TPL_INCLUDES})
    endif()
    if(TPL_COMPILE_DEFINITIONS)
      target_compile_definitions(${NAME} INTERFACE ${TPL_COMPILE_DEFINITIONS})
    endif()
    if(TPL_COMPILE_OPTIONS)
      target_compile_options(${NAME} INTERFACE ${TPL_COMPILE_OPTIONS})
    endif()
    if(TPL_LINK_OPTIONS)
      target_link_libraries(${NAME} INTERFACE ${TPL_LINK_OPTIONS})
    endif()
  else()
    add_library(${NAME} UNKNOWN IMPORTED)
    if(TPL_LIBRARY)
      set_target_properties(${NAME} PROPERTIES IMPORTED_LOCATION ${TPL_LIBRARY})
    endif()
    #Things have to go in quoted in case we have multiple list entries
    if(TPL_LINK_LIBRARIES)
      set_target_properties(${NAME} PROPERTIES INTERFACE_LINK_LIBRARIES "${TPL_LINK_LIBRARIES}")
    endif()
    if(TPL_INCLUDES)
      set_target_properties(${NAME} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${TPL_INCLUDES}")
    endif()
    if(TPL_COMPILE_DEFINITIONS)
      set_target_properties(${NAME} PROPERTIES INTERFACE_COMPILE_DEFINITIONS "${TPL_COMPILE_DEFINITIONS}")
    endif()
    if(TPL_COMPILE_OPTIONS)
      set_target_properties(${NAME} PROPERTIES INTERFACE_COMPILE_OPTIONS "${TPL_COMPILE_OPTIONS}")
    endif()
    if(TPL_LINK_OPTIONS)
      set_target_properties(${NAME} PROPERTIES INTERFACE_LINK_LIBRARIES "${TPL_LINK_OPTIONS}")
    endif()
  endif()
endmacro()

#
# @MACRO: KOKKOS_FIND_HEADER
#
# Function that finds a particular header. This searches custom paths
# or default system paths depending on options. In constrast to CMake
# default, custom paths are prioritized over system paths. The searched
# order is:
# 1. <NAME>_ROOT variable
# 2. <NAME>_ROOT environment variable
# 3. Kokkos_<NAME>_DIR variable
# 4. Locations in the PATHS option
# 5. Default system paths, if allowed.
#
# Default system paths are allowed if none of options (1)-(4) are specified
# or if default paths are specifically allowed via ALLOW_SYSTEM_PATH_FALLBACK
#
# Usage::
#
#   KOKKOS_FIND_HEADER(
#     <VAR_NAME>
#     <HEADER>
#     <TPL_NAME>
#    [ALLOW_SYSTEM_PATH_FALLBACK]
#    [PATHS path1 [path2 ...]]
#   )
#
#   ``<VAR_NAME>``
#
#   The variable to define with the success or failure of the find
#
#   ``<HEADER>``
#
#   The name of the header to find
#
#   ``<TPL_NAME>``
#
#   The name of the TPL the header corresponds to
#
#   ``[ALLOW_SYSTEM_PATH_FALLBACK]``
#
#   If custom paths are given and the header is not found
#   should we be allowed to search default system paths
#   or error out if not found in given paths
#
#   ``[PATHS path1 [path2 ...]]``
#
#   Custom paths to search for the header
#
macro(kokkos_find_header VAR_NAME HEADER TPL_NAME)
  cmake_parse_arguments(TPL "ALLOW_SYSTEM_PATH_FALLBACK" "" "PATHS" ${ARGN})

  set(${VAR_NAME} "${VARNAME}-NOTFOUND")
  set(HAVE_CUSTOM_PATHS FALSE)

  if(DEFINED ${TPL_NAME}_ROOT
     OR DEFINED ENV{${TPL_NAME}_ROOT}
     OR DEFINED KOKKOS_${TPL_NAME}_DIR
     OR TPL_PATHS
  )
    find_path(
      ${VAR_NAME} ${HEADER}
      PATHS ${${TPL_NAME}_ROOT} $ENV{${TPL_NAME}_ROOT} ${KOKKOS_${TPL_NAME}_DIR} ${TPL_PATHS}
      PATH_SUFFIXES include
      NO_DEFAULT_PATH
    )
    set(HAVE_CUSTOM_PATHS TRUE)
  endif()

  if(NOT HAVE_CUSTOM_PATHS OR TPL_ALLOW_SYSTEM_PATH_FALLBACK)
    #No-op if ${VAR_NAME} set by previous call
    find_path(${VAR_NAME} ${HEADER})
  endif()

endmacro()

#
# @MACRO: KOKKOS_FIND_LIBRARY
#
# Function that find a particular library. This searches custom paths
# or default system paths depending on options. In constrast to CMake
# default, custom paths are prioritized over system paths. The search
# order is:
# 1. <NAME>_ROOT variable
# 2. <NAME>_ROOT environment variable
# 3. Kokkos_<NAME>_DIR variable
# 4. Locations in the PATHS option
# 5. Default system paths, if allowed.
#
# Default system paths are allowed if none of options (1)-(3) are specified
# or if default paths are specifically allowed via ALLOW_SYSTEM_PATH_FALLBACK
#
# Usage::
#
#   KOKKOS_FIND_LIBRARY(
#     <VAR_NAME>
#     <HEADER>
#     <TPL_NAME>
#    [ALLOW_SYSTEM_PATH_FALLBACK]
#    [PATHS path1 [path2 ...]]
#    [SUFFIXES suffix1 [suffix2 ...]]
#   )
#
#   ``<VAR_NAME>``
#
#   The variable to define with the success or failure of the find
#
#   ``<LIBRARY>``
#
#   The name of the library to find (NOT prefixed with -l)
#
#   ``<TPL_NAME>``
#
#   The name of the TPL the library corresponds to
#
#   ``ALLOW_SYSTEM_PATH_FALLBACK``
#
#   If custom paths are given and the library is not found
#   should we be allowed to search default system paths
#   or error out if not found in given paths
#
#   ``PATHS``
#
#   Custom paths to search for the library
#
#   ``SUFFIXES``
#
#   Suffixes appended to PATHS when attempting to locate
#   the library. Defaults to {lib, lib64}.
#
macro(kokkos_find_library VAR_NAME LIB TPL_NAME)
  cmake_parse_arguments(TPL "ALLOW_SYSTEM_PATH_FALLBACK" "" "PATHS;SUFFIXES" ${ARGN})

  if(NOT TPL_SUFFIXES)
    set(TPL_SUFFIXES lib lib64)
  endif()

  set(${VAR_NAME} "${VARNAME}-NOTFOUND")
  set(HAVE_CUSTOM_PATHS FALSE)

  if(DEFINED ${TPL_NAME}_ROOT
     OR DEFINED ENV{${TPL_NAME}_ROOT}
     OR DEFINED KOKKOS_${TPL_NAME}_DIR
     OR TPL_PATHS
  )
    find_library(
      ${VAR_NAME} ${LIB}
      PATHS ${${TPL_NAME}_ROOT} $ENV{${TPL_NAME}_ROOT} ${KOKKOS_${TPL_NAME}_DIR} ${TPL_PATHS}
      PATH_SUFFIXES ${TPL_SUFFIXES}
      NO_DEFAULT_PATH
    )
    set(HAVE_CUSTOM_PATHS TRUE)
  endif()

  if(NOT HAVE_CUSTOM_PATHS OR TPL_ALLOW_SYSTEM_PATH_FALLBACK)
    #No-op if ${VAR_NAME} set by previous call
    find_library(${VAR_NAME} ${LIB} PATH_SUFFIXES ${TPL_SUFFIXES})
  endif()

endmacro()

#
# @MACRO: KOKKOS_FIND_IMPORTED
#
# Function that finds all libraries and headers needed for the tpl
# and creates an imported target encapsulating all the flags and libraries
#
# Usage::
#
#   KOKKOS_FIND_IMPORTED(
#     <NAME>
#     INTERFACE
#     ALLOW_SYSTEM_PATH_FALLBACK
#     MODULE_NAME <name>
#     IMPORTED_NAME <name>
#     LIBRARY <name>
#     LIBRARIES <name1> <name2> ...
#     LIBRARY_PATHS <path1> <path2> ...
#     LIBRARY_SUFFIXES <suffix1> <suffix2> ...
#     HEADER <name>
#     HEADERS <name1> <name2> ...
#     HEADER_PATHS <path1> <path2> ...
#   )
#
#   ``INTERFACE``
#
#     If specified, this TPL will build an INTERFACE library rather than an
#     IMPORTED target
#
#   ``ALLOW_SYSTEM_PATH_FALLBACK``
#
#     If custom paths are given and the library is not found
#     should we be allowed to search default system paths
#     or error out if not found in given paths.
#
#   ``MODULE_NAME <name>``
#
#     If specified, the name of the enclosing module passed to
#     FIND_PACKAGE(<MODULE_NAME>). Defaults to TPL${NAME} if not
#     given.
#
#   ``IMPORTED_NAME <name>``
#
#     If specified, this gives the name of the target to build.
#     Defaults to Kokkos::<NAME>
#
#   ``LIBRARY <name>``
#
#     If specified, this gives the name of the library to look for.
#     The full path for the library found will be used as IMPORTED_LOCATION
#     for the target created. Thus, this cannot be used for interface libraries.
#
#   ``LIBRARIES <name1> <name2> ...``
#
#     If specified, this gives a list of libraries to find for the package.
#     As opposed to the LIBRARY argument, this can be used with interface
#     libraries. In that case, we directly use the names provided here
#     for linking when creating the new target.
#
#   ``LIBRARY_PATHS <path1> <path2> ...``
#
#     If specified, this gives a list of paths to search for the library.
#     If not given, <NAME>_ROOT will be searched.
#
#   ``LIBRARY_SUFFIXES <suffix1> <suffix2> ...``
#
#     Suffixes appended to LIBRARY_PATHS when attempting to locate
#     libraries. If not given, defaults to {lib, lib64}.
#
#   ``HEADER <name>``
#
#     If specified, this gives the name of a header to to look for
#
#   ``HEADERS <name1> <name2> ...``
#
#     If specified, this gives a list of headers to find for the package
#
#   ``HEADER_PATHS <path1> <path2> ...``
#
#     If specified, this gives a list of paths to search for the headers
#     If not given, <NAME>_ROOT/include and <NAME>_ROOT/include will be searched.
#
macro(kokkos_find_imported NAME)
  cmake_parse_arguments(
    TPL "INTERFACE;ALLOW_SYSTEM_PATH_FALLBACK" "IMPORTED_NAME;MODULE_NAME;LIBRARY;HEADER"
    "LIBRARIES;LIBRARY_PATHS;LIBRARY_SUFFIXES;HEADERS;HEADER_PATHS" ${ARGN}
  )

  if(NOT TPL_MODULE_NAME)
    set(TPL_MODULE_NAME TPL${NAME})
  endif()

  if(TPL_ALLOW_SYSTEM_PATH_FALLBACK)
    set(ALLOW_PATH_FALLBACK_OPT ALLOW_SYSTEM_PATH_FALLBACK)
  else()
    set(ALLOW_PATH_FALLBACK_OPT)
  endif()

  if(NOT TPL_IMPORTED_NAME)
    if(TPL_INTERFACE)
      set(TPL_IMPORTED_NAME ${NAME})
    else()
      set(TPL_IMPORTED_NAME Kokkos::${NAME})
    endif()
  endif()

  if(NOT TPL_LIBRARY_SUFFIXES)
    set(TPL_LIBRARY_SUFFIXES lib)
    if(KOKKOS_IMPL_32BIT)
      list(APPEND TPL_LIBRARY_SUFFIXES lib32)
    else()
      list(APPEND TPL_LIBRARY_SUFFIXES lib64)
    endif()
  endif()

  set(${NAME}_INCLUDE_DIRS)
  if(TPL_HEADER)
    kokkos_find_header(${NAME}_INCLUDE_DIRS ${TPL_HEADER} ${NAME} ${ALLOW_PATH_FALLBACK_OPT} PATHS ${TPL_HEADER_PATHS})
  endif()

  foreach(HEADER ${TPL_HEADERS})
    kokkos_find_header(HEADER_FIND_TEMP ${HEADER} ${NAME} ${ALLOW_PATH_FALLBACK_OPT} PATHS ${TPL_HEADER_PATHS})
    if(HEADER_FIND_TEMP)
      list(APPEND ${NAME}_INCLUDE_DIRS ${HEADER_FIND_TEMP})
    endif()
  endforeach()

  set(${NAME}_LIBRARY)
  if(TPL_LIBRARY)
    kokkos_find_library(
      ${NAME}_LIBRARY
      ${TPL_LIBRARY}
      ${NAME}
      ${ALLOW_PATH_FALLBACK_OPT}
      PATHS
      ${TPL_LIBRARY_PATHS}
      SUFFIXES
      ${TPL_LIBRARY_SUFFIXES}
    )
  endif()

  set(${NAME}_FOUND_LIBRARIES)
  foreach(LIB ${TPL_LIBRARIES})
    kokkos_find_library(
      ${LIB}_LOCATION
      ${LIB}
      ${NAME}
      ${ALLOW_PATH_FALLBACK_OPT}
      PATHS
      ${TPL_LIBRARY_PATHS}
      SUFFIXES
      ${TPL_LIBRARY_SUFFIXES}
    )
    if(${LIB}_LOCATION)
      list(APPEND ${NAME}_FOUND_LIBRARIES ${${LIB}_LOCATION})
    else()
      set(${NAME}_FOUND_LIBRARIES ${${LIB}_LOCATION})
      break()
    endif()
  endforeach()

  include(FindPackageHandleStandardArgs)
  #Collect all the variables we need to be valid for
  #find_package to have succeeded
  set(TPL_VARS_NEEDED)
  if(TPL_LIBRARY)
    list(APPEND TPL_VARS_NEEDED ${NAME}_LIBRARY)
  endif()
  if(TPL_HEADER)
    list(APPEND TPL_VARS_NEEDED ${NAME}_INCLUDE_DIRS)
  endif()
  if(TPL_LIBRARIES)
    list(APPEND TPL_VARS_NEEDED ${NAME}_FOUND_LIBRARIES)
  endif()
  find_package_handle_standard_args(${TPL_MODULE_NAME} REQUIRED_VARS ${TPL_VARS_NEEDED})

  mark_as_advanced(${NAME}_INCLUDE_DIRS ${NAME}_FOUND_LIBRARIES ${NAME}_LIBRARY)

  #this is so much fun on a Cray system
  #/usr/include should never be added as a -isystem include
  #this freaks out the compiler include search order
  if(KOKKOS_IS_CRAYPE)
    list(REMOVE_ITEM ${NAME}_INCLUDE_DIRS "/usr/include")
  endif()

  if(${TPL_MODULE_NAME}_FOUND)
    set(IMPORT_TYPE)
    if(TPL_INTERFACE)
      set(IMPORT_TYPE "INTERFACE")
      set(${NAME}_FOUND_LIBRARIES ${TPL_LIBRARIES})
    endif()
    kokkos_create_imported_tpl(
      ${TPL_IMPORTED_NAME}
      ${IMPORT_TYPE}
      INCLUDES
      "${${NAME}_INCLUDE_DIRS}"
      LIBRARY
      "${${NAME}_LIBRARY}"
      LINK_LIBRARIES
      "${${NAME}_FOUND_LIBRARIES}"
    )
  endif()
endmacro(kokkos_find_imported)

#
# @MACRO: KOKKOS_LINK_TPL()
#
# Function that checks if a third-party library (TPL) has been enabled and
# calls target_link_libraries on the given target
#
# Usage::
#
#   KOKKOS_LINK_TPL(
#     <TARGET>
#     PUBLIC
#     PRIVATE
#     INTERFACE
#     IMPORTED_NAME  <name>
#     <TPL_NAME>
#
#   Checks if Kokkos_ENABLE_<TPL_NAME>=ON and if so links the library
#
#   ``PUBLIC/PRIVATE/INTERFACE``
#
#     Specifies the linkage mode. One of these arguments should be given.
#     This will then invoke target_link_libraries(<TARGET> PUBLIC/PRIVATE/INTERFACE <TPL_NAME>)
#
#   ``IMPORTED_NAME <name>``
#
#     If specified, this gives the exact name of the target to link against
#     target_link_libraries(<TARGET> <IMPORTED_NAME>)
#
function(kokkos_link_tpl TARGET)
  cmake_parse_arguments(TPL "PUBLIC;PRIVATE;INTERFACE" "IMPORTED_NAME" "" ${ARGN})
  #the name of the TPL
  set(TPL ${TPL_UNPARSED_ARGUMENTS})
  if(NOT TPL_IMPORTED_NAME)
    set(TPL_IMPORTED_NAME Kokkos::${TPL})
  endif()
  if(KOKKOS_ENABLE_${TPL})
    if(TPL_PUBLIC)
      target_link_libraries(${TARGET} PUBLIC ${TPL_IMPORTED_NAME})
    elseif(TPL_PRIVATE)
      target_link_libraries(${TARGET} PRIVATE ${TPL_IMPORTED_NAME})
    elseif(TPL_INTERFACE)
      target_link_libraries(${TARGET} INTERFACE ${TPL_IMPORTED_NAME})
    else()
      target_link_libraries(${TARGET} ${TPL_IMPORTED_NAME})
    endif()
  endif()
endfunction()

function(COMPILER_SPECIFIC_OPTIONS_HELPER)
  set(COMPILERS
      NVIDIA
      NVHPC
      DEFAULT
      Cray
      Clang
      AppleClang
      IntelLLVM
      GNU
      HIPCC
      Fujitsu
      MSVC
      CrayClang
  )
  cmake_parse_arguments(
    PARSE "LINK_OPTIONS;COMPILE_OPTIONS;COMPILE_DEFINITIONS;LINK_LIBRARIES" "COMPILER_ID" "${COMPILERS}" ${ARGN}
  )
  if(PARSE_UNPARSED_ARGUMENTS)
    message(
      SEND_ERROR "'${PARSE_UNPARSED_ARGUMENTS}' argument(s) not recognized when providing compiler specific options"
    )
  endif()

  if(PARSE_COMPILER_ID)
    set(COMPILER ${${PARSE_COMPILER_ID}})
  else()
    set(COMPILER ${KOKKOS_CXX_COMPILER_ID})
  endif()

  set(COMPILER_SPECIFIC_FLAGS_TMP ${PARSE_DEFAULT})
  foreach(COMP ${COMPILERS})
    if(COMPILER STREQUAL "${COMP}")
      if(PARSE_${COMPILER})
        if("${PARSE_${COMPILER}}" STREQUAL "NO-VALUE-SPECIFIED")
          set(COMPILER_SPECIFIC_FLAGS_TMP "")
        else()
          set(COMPILER_SPECIFIC_FLAGS_TMP ${PARSE_${COMPILER}})
        endif()
      endif()
    endif()
  endforeach()

  if(PARSE_COMPILE_OPTIONS)
    # The funky logic here is for future handling of argument deduplication
    # If we naively pass multiple -Xcompiler flags to target_compile_options
    # -Xcompiler will get deduplicated and break the build
    if("-Xcompiler" IN_LIST COMPILER_SPECIFIC_FLAGS_TMP)
      list(REMOVE_ITEM COMPILER_SPECIFIC_FLAGS_TMP "-Xcompiler")
      global_append(KOKKOS_XCOMPILER_OPTIONS ${COMPILER_SPECIFIC_FLAGS_TMP})
    else()
      global_append(KOKKOS_COMPILE_OPTIONS ${COMPILER_SPECIFIC_FLAGS_TMP})
    endif()
  endif()

  if(PARSE_LINK_OPTIONS)
    global_append(KOKKOS_LINK_OPTIONS ${COMPILER_SPECIFIC_FLAGS_TMP})
  endif()

  if(PARSE_COMPILE_DEFINITIONS)
    global_append(KOKKOS_COMPILE_DEFINITIONS ${COMPILER_SPECIFIC_FLAGS_TMP})
  endif()

  if(PARSE_LINK_LIBRARIES)
    global_append(KOKKOS_LINK_LIBRARIES ${COMPILER_SPECIFIC_FLAGS_TMP})
  endif()
endfunction(COMPILER_SPECIFIC_OPTIONS_HELPER)

function(COMPILER_SPECIFIC_FLAGS)
  compiler_specific_options_helper(${ARGN} COMPILE_OPTIONS LINK_OPTIONS)
endfunction(COMPILER_SPECIFIC_FLAGS)

function(COMPILER_SPECIFIC_OPTIONS)
  compiler_specific_options_helper(${ARGN} COMPILE_OPTIONS)
endfunction(COMPILER_SPECIFIC_OPTIONS)

function(COMPILER_SPECIFIC_LINK_OPTIONS)
  compiler_specific_options_helper(${ARGN} LINK_OPTIONS)
endfunction(COMPILER_SPECIFIC_LINK_OPTIONS)

function(COMPILER_SPECIFIC_DEFS)
  compiler_specific_options_helper(${ARGN} COMPILE_DEFINITIONS)
endfunction(COMPILER_SPECIFIC_DEFS)

function(COMPILER_SPECIFIC_LIBS)
  compiler_specific_options_helper(${ARGN} LINK_LIBRARIES)
endfunction(COMPILER_SPECIFIC_LIBS)
# Given a list of the form
#  key1;value1;key2;value2,...
# Create a list of all keys in a variable named ${KEY_LIST_NAME}
# and set the value for each key in a variable ${VAR_PREFIX}key1,...
# kokkos_key_value_map(ARCH ALL_ARCHES key1;value1;key2;value2)
# would produce a list variable ALL_ARCHES=key1;key2
# and individual variables ARCHkey1=value1 and ARCHkey2=value2
macro(KOKKOS_KEY_VALUE_MAP VAR_PREFIX KEY_LIST_NAME)
  set(PARSE_KEY ON)
  set(${KEY_LIST_NAME})
  foreach(ENTRY ${ARGN})
    if(PARSE_KEY)
      set(CURRENT_KEY ${ENTRY})
      set(PARSE_KEY OFF)
      list(APPEND ${KEY_LIST_NAME} ${CURRENT_KEY})
    else()
      set(${VAR_PREFIX}${CURRENT_KEY} ${ENTRY})
      set(PARSE_KEY ON)
    endif()
  endforeach()
endmacro()

function(KOKKOS_CHECK_DEPRECATED_OPTIONS)
  kokkos_key_value_map(DEPRECATED_MSG_ DEPRECATED_LIST ${ARGN})
  foreach(OPTION_SUFFIX ${DEPRECATED_LIST})
    set(OPTION_NAME Kokkos_${OPTION_SUFFIX})
    set(OPTION_MESSAGE ${DEPRECATED_MSG_${OPTION_SUFFIX}})
    if(DEFINED ${OPTION_NAME}) # This variable has been given by the user as on or off
      message(SEND_ERROR "Removed option ${OPTION_NAME} has been given with value ${${OPTION_NAME}}. ${OPT_MESSAGE}")
    endif()
  endforeach()
endfunction()

# this function checks whether the current CXX compiler supports building CUDA
function(kokkos_cxx_compiler_cuda_test _VAR)
  # don't run this test every time
  if(DEFINED ${_VAR})
    return()
  endif()

  file(
    WRITE ${PROJECT_BINARY_DIR}/compile_tests/compiles_cuda.cpp
    "
#include <cuda.h>
#include <cstdlib>

__global__
void kernel(int sz, double* data)
{
    auto _beg = blockIdx.x * blockDim.x + threadIdx.x;
    for(int i = _beg; i < sz; ++i)
        data[i] += static_cast<double>(i);
}

int main()
{
    double* data = nullptr;
    int blocks = 64;
    int grids = 64;
    auto ret = cudaMalloc(&data, blocks * grids * sizeof(double));
    if(ret != cudaSuccess)
        return EXIT_FAILURE;
    kernel<<<grids, blocks>>>(blocks * grids, data);
    cudaDeviceSynchronize();
    return EXIT_SUCCESS;
}
"
  )

  try_compile(_RET ${PROJECT_BINARY_DIR}/compile_tests SOURCES ${PROJECT_BINARY_DIR}/compile_tests/compiles_cuda.cpp)

  set(${_VAR} ${_RET} CACHE STRING "CXX compiler supports building CUDA")
endfunction()

# this function is provided to easily select which files use nvcc_wrapper:
#
#       COMPILER    --> do check for compiler
#       LINKER      --> do check for linker
#       LANGUAGE    --> the language for which to check (required)
#       FLAGS       --> flags to check
#
function(kokkos_check_flags)
  cmake_parse_arguments(INP "COMPILER;LINKER" "LANGUAGE" "FLAGS;LINKER_FLAGS" ${ARGN})

  # do nothing if no flags are given
  if(NOT INP_FLAGS)
    return()
  endif()

  if(NOT INP_LANGUAGE)
    message(FATAL_ERROR "'kokkos_check_flags' requires LANGUAGE keyword")
  endif()

  #check_compiler/linker_flag requires a whitespace separated list
  string(REPLACE ";" " " WHITESPACE_FLAGS "${INP_FLAGS}")
  # icpx adds "-device ..." options that need quotes (which CMake removes). We need to add them here again.
  string(REGEX REPLACE "(-device [A-Za-z0-9_\\\\.]*)" "\"\\1\"" QUOTED_FLAGS "${WHITESPACE_FLAGS}")

  if(INP_COMPILER)
    include(CheckCompilerFlag)
    #delete cache so we always do the check
    unset(KOKKOS_COMPILE_OPTIONS_CHECK CACHE)
    if(INP_LINKER_FLAGS)
      # icpx adds "-device ..." options that need quotes (which CMake removes). We need to add them here again.
      string(REGEX REPLACE "(-device [A-Za-z0-9_\\\\.]*)" "\"\\1\"" QUOTED_LINKER_FLAGS "${INP_LINKER_FLAGS}")
      set(CMAKE_REQUIRED_LINK_OPTIONS "${QUOTED_LINKER_FLAGS}")
    endif()
    check_compiler_flag(${INP_LANGUAGE} "${QUOTED_FLAGS}" KOKKOS_COMPILE_OPTIONS_CHECK)
    if(NOT KOKKOS_COMPILE_OPTIONS_CHECK)
      message(
        FATAL_ERROR
          "The compiler for ${KOKKOS_COMPILE_LANGUAGE} can not consume flag(s) ${QUOTED_FLAGS} in combination with the CMAKE_${KOKKOS_COMPILE_LANGUAGE}_FLAGS=${CMAKE_${KOKKOS_COMPILE_LANGUAGE}_FLAGS}. Please check the given configuration."
      )
    endif()
  endif()

  if(INP_LINKER)
    include(CheckLinkerFlag)
    # temporarily set language flags to nothing ... the linker often cannot handle these which leads to false errors
    set(CMAKE_${INP_LANGUAGE}_FLAGS "")
    #delete cache so we always do the check
    unset(KOKKOS_LINK_OPTIONS_CHECK CACHE)
    check_linker_flag(${INP_LANGUAGE} "${QUOTED_FLAGS}" KOKKOS_LINK_OPTIONS_CHECK)
    if(NOT KOKKOS_LINK_OPTIONS_CHECK)
      message(
        FATAL_ERROR
          "The linker for ${KOKKOS_COMPILE_LANGUAGE} can not consume flag(s) ${QUOTED_FLAGS}. Please check the given configuration."
      )
    endif()
  endif()
endfunction()

# this function is provided to easily select which files use nvcc_wrapper:
#
#       GLOBAL      --> all files
#       TARGET      --> all files in a target
#       SOURCE      --> specific source files
#       DIRECTORY   --> all files in directory
#       PROJECT     --> all files/targets in a project/subproject
#
# NOTE: this is VERY DIFFERENT than the version in KokkosConfigCommon.cmake.in.
# This version explicitly uses nvcc_wrapper.
#
function(kokkos_compilation)
  # check whether the compiler already supports building CUDA
  kokkos_cxx_compiler_cuda_test(Kokkos_CXX_COMPILER_COMPILES_CUDA)
  # if CUDA compile test has already been performed, just return
  if(Kokkos_CXX_COMPILER_COMPILES_CUDA)
    return()
  endif()

  cmake_parse_arguments(COMP "GLOBAL;PROJECT" "" "DIRECTORY;TARGET;SOURCE" ${ARGN})

  # find kokkos_launch_compiler
  find_program(
    Kokkos_COMPILE_LAUNCHER
    NAMES kokkos_launch_compiler
    HINTS ${PROJECT_SOURCE_DIR}
    PATHS ${PROJECT_SOURCE_DIR}
    PATH_SUFFIXES bin
  )

  if(NOT Kokkos_COMPILE_LAUNCHER)
    message(
      FATAL_ERROR
        "Kokkos could not find 'kokkos_launch_compiler'. Please set '-DKokkos_COMPILE_LAUNCHER=/path/to/launcher'"
    )
  endif()

  # find nvcc_wrapper
  find_program(
    Kokkos_NVCC_WRAPPER
    NAMES nvcc_wrapper
    HINTS ${PROJECT_SOURCE_DIR}
    PATHS ${PROJECT_SOURCE_DIR}
    PATH_SUFFIXES bin
  )

  if(NOT Kokkos_COMPILE_LAUNCHER)
    message(
      FATAL_ERROR "Kokkos could not find 'nvcc_wrapper'. Please set '-DKokkos_COMPILE_LAUNCHER=/path/to/nvcc_wrapper'"
    )
  endif()

  if(COMP_GLOBAL)
    # if global, don't bother setting others
    set_property(
      GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${Kokkos_COMPILE_LAUNCHER} ${Kokkos_NVCC_WRAPPER} ${CMAKE_CXX_COMPILER}"
    )
    set_property(
      GLOBAL PROPERTY RULE_LAUNCH_LINK "${Kokkos_COMPILE_LAUNCHER} ${Kokkos_NVCC_WRAPPER} ${CMAKE_CXX_COMPILER}"
    )
  else()
    foreach(_TYPE PROJECT DIRECTORY TARGET SOURCE)
      # make project/subproject scoping easy, e.g. KokkosCompilation(PROJECT) after project(...)
      if("${_TYPE}" STREQUAL "PROJECT" AND COMP_${_TYPE})
        list(APPEND COMP_DIRECTORY ${PROJECT_SOURCE_DIR})
        unset(COMP_${_TYPE})
      endif()
      # set the properties if defined
      if(COMP_${_TYPE})
        # MESSAGE(STATUS "Using nvcc_wrapper :: ${_TYPE} :: ${COMP_${_TYPE}}")
        set_property(
          ${_TYPE} ${COMP_${_TYPE}} PROPERTY RULE_LAUNCH_COMPILE
                                             "${Kokkos_COMPILE_LAUNCHER} ${Kokkos_NVCC_WRAPPER} ${CMAKE_CXX_COMPILER}"
        )
        set_property(
          ${_TYPE} ${COMP_${_TYPE}} PROPERTY RULE_LAUNCH_LINK
                                             "${Kokkos_COMPILE_LAUNCHER} ${Kokkos_NVCC_WRAPPER} ${CMAKE_CXX_COMPILER}"
        )
      endif()
    endforeach()
  endif()
endfunction()
## KOKKOS_CONFIG_HEADER - parse the data list which is a list of backend names
##                        and create output config header file...used for
##                        creating dynamic include files based on enabled backends
##
##                        SRC_FILE is input file
##                        TARGET_FILE output file
##                        HEADER_GUARD TEXT used with include header guard
##                        HEADER_PREFIX prefix used with include (i.e. fwd, decl, setup)
##                        DATA_LIST list of backends to include in generated file
function(KOKKOS_CONFIG_HEADER SRC_FILE TARGET_FILE HEADER_GUARD HEADER_PREFIX DATA_LIST)
  set(HEADER_GUARD_TAG "${HEADER_GUARD}_HPP_")
  configure_file(cmake/${SRC_FILE} ${PROJECT_BINARY_DIR}/temp/${TARGET_FILE}.work COPYONLY)
  foreach(BACKEND_NAME ${DATA_LIST})
    set(INCLUDE_NEXT_FILE "#include <${HEADER_PREFIX}_${BACKEND_NAME}.hpp>
\@INCLUDE_NEXT_FILE\@"
    )
    configure_file(${PROJECT_BINARY_DIR}/temp/${TARGET_FILE}.work ${PROJECT_BINARY_DIR}/temp/${TARGET_FILE}.work @ONLY)
  endforeach()
  set(INCLUDE_NEXT_FILE "")
  configure_file(${PROJECT_BINARY_DIR}/temp/${TARGET_FILE}.work ${TARGET_FILE} @ONLY)
endfunction()
