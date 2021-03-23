################################### FUNCTIONS ##################################
# List of functions
#   kokkos_option

# Validate options are given with correct case and define an internal
# upper-case version for use within

#
#
# @FUNCTION: kokkos_deprecated_list
#
# Function that checks if a deprecated list option like Kokkos_ARCH was given.
# This prints an error and prevents configure from completing.
# It attempts to print a helpful message about updating the options for the new CMake.
# Kokkos_${SUFFIX} is the name of the option (like Kokkos_ARCH) being checked.
# Kokkos_${PREFIX}_X is the name of new option to be defined from a list X,Y,Z,...
FUNCTION(kokkos_deprecated_list SUFFIX PREFIX)
  SET(CAMEL_NAME Kokkos_${SUFFIX})
  STRING(TOUPPER ${CAMEL_NAME} UC_NAME)

  #I don't love doing it this way but better to be safe
  FOREACH(opt ${KOKKOS_GIVEN_VARIABLES})
    STRING(TOUPPER ${opt} OPT_UC)
    IF ("${OPT_UC}" STREQUAL "${UC_NAME}")
      STRING(REPLACE "," ";" optlist "${${opt}}")
      SET(ERROR_MSG "Given deprecated option list ${opt}. This must now be given as separate -D options, which assuming you spelled options correctly would be:")
      FOREACH(entry ${optlist})
        STRING(TOUPPER ${entry} ENTRY_UC)
        STRING(APPEND ERROR_MSG "\n  -DKokkos_${PREFIX}_${ENTRY_UC}=ON")
      ENDFOREACH()
      STRING(APPEND ERROR_MSG "\nRemove CMakeCache.txt and re-run. For a list of valid options, refer to BUILD.md or even look at CMakeCache.txt (before deleting it).")
      IF (KOKKOS_HAS_TRILINOS)
        MESSAGE(WARNING ${ERROR_MSG})
        FOREACH(entry ${optlist})
          STRING(TOUPPER ${entry} ENTRY_UC)
          SET(${CAMEL_NAME}_${ENTRY_UC} ON CACHE BOOL "Deprecated Trilinos translation")
        ENDFOREACH()
        UNSET(${opt} CACHE)
      ELSE()
        MESSAGE(SEND_ERROR ${ERROR_MSG})
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDFUNCTION()

FUNCTION(kokkos_option CAMEL_SUFFIX DEFAULT TYPE DOCSTRING)
  SET(CAMEL_NAME Kokkos_${CAMEL_SUFFIX})
  STRING(TOUPPER ${CAMEL_NAME} UC_NAME)

  LIST(APPEND KOKKOS_OPTION_KEYS ${CAMEL_SUFFIX})
  SET(KOKKOS_OPTION_KEYS ${KOKKOS_OPTION_KEYS} PARENT_SCOPE)
  LIST(APPEND KOKKOS_OPTION_VALUES "${DOCSTRING}")
  SET(KOKKOS_OPTION_VALUES ${KOKKOS_OPTION_VALUES} PARENT_SCOPE)
  LIST(APPEND KOKKOS_OPTION_TYPES ${TYPE})
  SET(KOKKOS_OPTION_TYPES ${KOKKOS_OPTION_TYPES} PARENT_SCOPE)

  # Make sure this appears in the cache with the appropriate DOCSTRING
  SET(${CAMEL_NAME} ${DEFAULT} CACHE ${TYPE} ${DOCSTRING})

  #I don't love doing it this way because it's N^2 in number options, but cest la vie
  FOREACH(opt ${KOKKOS_GIVEN_VARIABLES})
    STRING(TOUPPER ${opt} OPT_UC)
    IF ("${OPT_UC}" STREQUAL "${UC_NAME}")
      IF (NOT "${opt}" STREQUAL "${CAMEL_NAME}")
        IF (KOKKOS_HAS_TRILINOS)
          #Allow this for now if Trilinos... we need to bootstrap our way to integration
          MESSAGE(WARNING "Deprecated option ${opt} found - please change spelling to ${CAMEL_NAME}")
          SET(${CAMEL_NAME} "${${opt}}" CACHE ${TYPE} ${DOCSTRING} FORCE)
          UNSET(${opt} CACHE)
        ELSE()
          MESSAGE(FATAL_ERROR "Matching option found for ${CAMEL_NAME} with the wrong case ${opt}. Please delete your CMakeCache.txt and change option to -D${CAMEL_NAME}=${${opt}}. This is now enforced to avoid hard-to-debug CMake cache inconsistencies.")
        ENDIF()
      ENDIF()
    ENDIF()
  ENDFOREACH()

  #okay, great, we passed the validation test - use the default
  IF (DEFINED ${CAMEL_NAME})
    SET(${UC_NAME} ${${CAMEL_NAME}} PARENT_SCOPE)
  ELSE()
    SET(${UC_NAME} ${DEFAULT} PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

FUNCTION(kokkos_set_option CAMEL_SUFFIX VALUE)
  LIST(FIND KOKKOS_OPTION_KEYS ${CAMEL_SUFFIX} OPTION_INDEX)
  IF(OPTION_INDEX EQUAL -1)
    MESSAGE(FATAL_ERROR "Couldn't set value for Kokkos_${CAMEL_SUFFIX}")
  ENDIF()
  SET(CAMEL_NAME Kokkos_${CAMEL_SUFFIX})
  STRING(TOUPPER ${CAMEL_NAME} UC_NAME)

  LIST(GET KOKKOS_OPTION_VALUES ${OPTION_INDEX} DOCSTRING)
  LIST(GET KOKKOS_OPTION_TYPES ${OPTION_INDEX} TYPE)
  SET(${CAMEL_NAME} ${VALUE} CACHE ${TYPE} ${DOCSTRING} FORCE)
  MESSAGE(STATUS "Setting ${CAMEL_NAME}=${VALUE}")
  SET(${UC_NAME} ${VALUE} PARENT_SCOPE)
ENDFUNCTION()

FUNCTION(kokkos_append_config_line LINE)
  GLOBAL_APPEND(KOKKOS_TPL_EXPORTS "${LINE}")
ENDFUNCTION()

MACRO(kokkos_export_cmake_tpl NAME)
  #CMake TPLs are located with a call to find_package
  #find_package locates XConfig.cmake files through
  #X_DIR or X_ROOT variables set prior to calling find_package

  #If Kokkos was configured to find the TPL through a _DIR variable
  #make sure thar DIR variable is available to downstream packages
  IF (DEFINED ${NAME}_DIR)
    #The downstream project may override the TPL location that Kokkos used
    #Check if the downstream project chose its own TPL location
    #If not, make the Kokkos found location available
    KOKKOS_APPEND_CONFIG_LINE("IF(NOT DEFINED ${NAME}_DIR)")
    KOKKOS_APPEND_CONFIG_LINE("  SET(${NAME}_DIR  ${${NAME}_DIR})")
    KOKKOS_APPEND_CONFIG_LINE("ENDIF()")
  ENDIF()

  IF (DEFINED ${NAME}_ROOT)
    #The downstream project may override the TPL location that Kokkos used
    #Check if the downstream project chose its own TPL location
    #If not, make the Kokkos found location available
    KOKKOS_APPEND_CONFIG_LINE("IF(NOT DEFINED ${NAME}_ROOT)")
    KOKKOS_APPEND_CONFIG_LINE("  SET(${NAME}_ROOT  ${${NAME}_ROOT})")
    KOKKOS_APPEND_CONFIG_LINE("ENDIF()")
  ENDIF()
  KOKKOS_APPEND_CONFIG_LINE("FIND_DEPENDENCY(${NAME})")
ENDMACRO()

MACRO(kokkos_export_imported_tpl NAME)
  IF (NOT KOKKOS_HAS_TRILINOS)
    GET_TARGET_PROPERTY(LIB_IMPORTED ${NAME} IMPORTED)
    IF (NOT LIB_IMPORTED)
      # This is not an imported target
      # This an interface library that we created
      INSTALL(
        TARGETS ${NAME}
        EXPORT KokkosTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      )
    ELSE()
      #make sure this also gets "exported" in the config file
      KOKKOS_APPEND_CONFIG_LINE("IF(NOT TARGET ${NAME})")

      GET_TARGET_PROPERTY(LIB_TYPE ${NAME} TYPE)
      IF (${LIB_TYPE} STREQUAL "INTERFACE_LIBRARY")
        KOKKOS_APPEND_CONFIG_LINE("ADD_LIBRARY(${NAME} INTERFACE IMPORTED)")
        KOKKOS_APPEND_CONFIG_LINE("SET_TARGET_PROPERTIES(${NAME} PROPERTIES")
      ELSE()
        KOKKOS_APPEND_CONFIG_LINE("ADD_LIBRARY(${NAME} UNKNOWN IMPORTED)")
        KOKKOS_APPEND_CONFIG_LINE("SET_TARGET_PROPERTIES(${NAME} PROPERTIES")
        GET_TARGET_PROPERTY(TPL_LIBRARY ${NAME} IMPORTED_LOCATION)
        IF(TPL_LIBRARY)
          KOKKOS_APPEND_CONFIG_LINE("IMPORTED_LOCATION \"${TPL_LIBRARY}\"")
        ENDIF()
      ENDIF()

      GET_TARGET_PROPERTY(TPL_INCLUDES ${NAME} INTERFACE_INCLUDE_DIRECTORIES)
      IF(TPL_INCLUDES)
        KOKKOS_APPEND_CONFIG_LINE("INTERFACE_INCLUDE_DIRECTORIES \"${TPL_INCLUDES}\"")
      ENDIF()

      GET_TARGET_PROPERTY(TPL_COMPILE_OPTIONS ${NAME} INTERFACE_COMPILE_OPTIONS)
      IF(TPL_COMPILE_OPTIONS)
        KOKKOS_APPEND_CONFIG_LINE("INTERFACE_COMPILE_OPTIONS ${TPL_COMPILE_OPTIONS}")
      ENDIF()

      SET(TPL_LINK_OPTIONS)
      IF(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13.0")
        GET_TARGET_PROPERTY(TPL_LINK_OPTIONS ${NAME} INTERFACE_LINK_OPTIONS)
      ENDIF()
      IF(TPL_LINK_OPTIONS)
        KOKKOS_APPEND_CONFIG_LINE("INTERFACE_LINK_OPTIONS ${TPL_LINK_OPTIONS}")
      ENDIF()

      GET_TARGET_PROPERTY(TPL_LINK_LIBRARIES  ${NAME} INTERFACE_LINK_LIBRARIES)
      IF(TPL_LINK_LIBRARIES)
        KOKKOS_APPEND_CONFIG_LINE("INTERFACE_LINK_LIBRARIES \"${TPL_LINK_LIBRARIES}\"")
      ENDIF()
      KOKKOS_APPEND_CONFIG_LINE(")")
      KOKKOS_APPEND_CONFIG_LINE("ENDIF()")
    ENDIF()
  ENDIF()
ENDMACRO()


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
IF (KOKKOS_HAS_TRILINOS)
MACRO(kokkos_import_tpl NAME)
  #do nothing
ENDMACRO()
ELSE()
MACRO(kokkos_import_tpl NAME)
  CMAKE_PARSE_ARGUMENTS(TPL
   "NO_EXPORT;INTERFACE"
   ""
   ""
   ${ARGN})
  IF (TPL_INTERFACE)
    SET(TPL_IMPORTED_NAME ${NAME})
  ELSE()
    SET(TPL_IMPORTED_NAME Kokkos::${NAME})
  ENDIF()

  # Even though this policy gets set in the top-level CMakeLists.txt,
  # I have still been getting errors about ROOT variables being ignored
  # I'm not sure if this is a scope issue - but make sure
  # the policy is set before we do any find_package calls
  IF(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
    CMAKE_POLICY(SET CMP0074 NEW)
  ENDIF()

  IF (KOKKOS_ENABLE_${NAME})
    #Tack on a TPL here to make sure we avoid using anyone else's find
    FIND_PACKAGE(TPL${NAME} REQUIRED MODULE)
    IF(NOT TARGET ${TPL_IMPORTED_NAME})
      MESSAGE(FATAL_ERROR "Find module succeeded for ${NAME}, but did not produce valid target ${TPL_IMPORTED_NAME}")
    ENDIF()
    IF(NOT TPL_NO_EXPORT)
      KOKKOS_EXPORT_IMPORTED_TPL(${TPL_IMPORTED_NAME})
    ENDIF()
    LIST(APPEND KOKKOS_ENABLED_TPLS ${NAME})
  ENDIF()
ENDMACRO(kokkos_import_tpl)
ENDIF()

MACRO(kokkos_import_cmake_tpl MODULE_NAME)
  kokkos_import_tpl(${MODULE_NAME} ${ARGN} NO_EXPORT)
  CMAKE_PARSE_ARGUMENTS(TPL
   "NO_EXPORT"
   "OPTION_NAME"
   ""
   ${ARGN})

  IF (NOT TPL_OPTION_NAME)
    SET(TPL_OPTION_NAME ${MODULE_NAME})
  ENDIF()

  IF (NOT TPL_NO_EXPORT)
    KOKKOS_EXPORT_CMAKE_TPL(${MODULE_NAME})
  ENDIF()
ENDMACRO()

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
MACRO(kokkos_create_imported_tpl NAME)
  CMAKE_PARSE_ARGUMENTS(TPL
   "INTERFACE"
   "LIBRARY"
   "LINK_LIBRARIES;INCLUDES;COMPILE_OPTIONS;LINK_OPTIONS"
   ${ARGN})


  IF (KOKKOS_HAS_TRILINOS)
    #TODO: we need to set a bunch of cache variables here
  ELSEIF (TPL_INTERFACE)
    ADD_LIBRARY(${NAME} INTERFACE)
    #Give this an importy-looking name
    ADD_LIBRARY(Kokkos::${NAME} ALIAS ${NAME})
    IF (TPL_LIBRARY)
      MESSAGE(SEND_ERROR "TPL Interface library ${NAME} should not have an IMPORTED_LOCATION")
    ENDIF()
    #Things have to go in quoted in case we have multiple list entries
    IF(TPL_LINK_LIBRARIES)
      TARGET_LINK_LIBRARIES(${NAME} INTERFACE ${TPL_LINK_LIBRARIES})
    ENDIF()
    IF(TPL_INCLUDES)
      TARGET_INCLUDE_DIRECTORIES(${NAME} INTERFACE ${TPL_INCLUDES})
    ENDIF()
    IF(TPL_COMPILE_OPTIONS)
      TARGET_COMPILE_OPTIONS(${NAME} INTERFACE ${TPL_COMPILE_OPTIONS})
    ENDIF()
    IF(TPL_LINK_OPTIONS)
      TARGET_LINK_LIBRARIES(${NAME} INTERFACE ${TPL_LINK_OPTIONS})
    ENDIF()
  ELSE()
    ADD_LIBRARY(${NAME} UNKNOWN IMPORTED)
    IF(TPL_LIBRARY)
      SET_TARGET_PROPERTIES(${NAME} PROPERTIES
        IMPORTED_LOCATION ${TPL_LIBRARY})
    ENDIF()
    #Things have to go in quoted in case we have multiple list entries
    IF(TPL_LINK_LIBRARIES)
      SET_TARGET_PROPERTIES(${NAME} PROPERTIES
        INTERFACE_LINK_LIBRARIES "${TPL_LINK_LIBRARIES}")
    ENDIF()
    IF(TPL_INCLUDES)
      SET_TARGET_PROPERTIES(${NAME} PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${TPL_INCLUDES}")
    ENDIF()
    IF(TPL_COMPILE_OPTIONS)
      SET_TARGET_PROPERTIES(${NAME} PROPERTIES
        INTERFACE_COMPILE_OPTIONS "${TPL_COMPILE_OPTIONS}")
    ENDIF()
    IF(TPL_LINK_OPTIONS)
      SET_TARGET_PROPERTIES(${NAME} PROPERTIES
        INTERFACE_LINK_LIBRARIES "${TPL_LINK_OPTIONS}")
    ENDIF()
  ENDIF()
ENDMACRO()

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
MACRO(kokkos_find_header VAR_NAME HEADER TPL_NAME)
  CMAKE_PARSE_ARGUMENTS(TPL
   "ALLOW_SYSTEM_PATH_FALLBACK"
   ""
   "PATHS"
   ${ARGN})

  SET(${VAR_NAME} "${VARNAME}-NOTFOUND")
  SET(HAVE_CUSTOM_PATHS FALSE)

  IF(DEFINED ${TPL_NAME}_ROOT OR
     DEFINED ENV{${TPL_NAME}_ROOT} OR
     DEFINED KOKKOS_${TPL_NAME}_DIR OR
     TPL_PATHS)
    FIND_PATH(${VAR_NAME} ${HEADER}
      PATHS
        ${${TPL_NAME}_ROOT}
        $ENV{${TPL_NAME}_ROOT}
        ${KOKKOS_${TPL_NAME}_DIR}
        ${TPL_PATHS}
      PATH_SUFFIXES include
      NO_DEFAULT_PATH)
    SET(HAVE_CUSTOM_PATHS TRUE)
  ENDIF()

  IF(NOT HAVE_CUSTOM_PATHS OR TPL_ALLOW_SYSTEM_PATH_FALLBACK)
    #No-op if ${VAR_NAME} set by previous call
    FIND_PATH(${VAR_NAME} ${HEADER})
  ENDIF()

ENDMACRO()

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
MACRO(kokkos_find_library VAR_NAME LIB TPL_NAME)
  CMAKE_PARSE_ARGUMENTS(TPL
   "ALLOW_SYSTEM_PATH_FALLBACK"
   ""
   "PATHS;SUFFIXES"
   ${ARGN})

  IF(NOT TPL_SUFFIXES)
    SET(TPL_SUFFIXES lib lib64)
  ENDIF()

  SET(${VAR_NAME} "${VARNAME}-NOTFOUND")
  SET(HAVE_CUSTOM_PATHS FALSE)

  IF(DEFINED ${TPL_NAME}_ROOT OR
     DEFINED ENV{${TPL_NAME}_ROOT} OR
     DEFINED KOKKOS_${TPL_NAME}_DIR OR
     TPL_PATHS)
    FIND_LIBRARY(${VAR_NAME} ${LIB}
      PATHS
        ${${TPL_NAME}_ROOT}
        $ENV{${TPL_NAME}_ROOT}
        ${KOKKOS_${TPL_NAME}_DIR}
        ${TPL_PATHS}
      PATH_SUFFIXES
        ${TPL_SUFFIXES}
      NO_DEFAULT_PATH)
    SET(HAVE_CUSTOM_PATHS TRUE)
  ENDIF()

  IF(NOT HAVE_CUSTOM_PATHS OR TPL_ALLOW_SYSTEM_PATH_FALLBACK)
    #No-op if ${VAR_NAME} set by previous call
    FIND_LIBRARY(${VAR_NAME} ${LIB} PATH_SUFFIXES ${TPL_SUFFIXES})
  ENDIF()

ENDMACRO()

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
#     If specified, this gives the name of the library to look for
#
#   ``LIBRARIES <name1> <name2> ...``
#
#     If specified, this gives a list of libraries to find for the package
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
MACRO(kokkos_find_imported NAME)
  CMAKE_PARSE_ARGUMENTS(TPL
   "INTERFACE;ALLOW_SYSTEM_PATH_FALLBACK"
   "IMPORTED_NAME;MODULE_NAME;LIBRARY;HEADER"
   "LIBRARIES;LIBRARY_PATHS;LIBRARY_SUFFIXES;HEADERS;HEADER_PATHS"
   ${ARGN})

  IF(NOT TPL_MODULE_NAME)
    SET(TPL_MODULE_NAME TPL${NAME})
  ENDIF()

  IF (TPL_ALLOW_SYSTEM_PATH_FALLBACK)
    SET(ALLOW_PATH_FALLBACK_OPT ALLOW_SYSTEM_PATH_FALLBACK)
  ELSE()
    SET(ALLOW_PATH_FALLBACK_OPT)
  ENDIF()

  IF (NOT TPL_IMPORTED_NAME)
    IF (TPL_INTERFACE)
      SET(TPL_IMPORTED_NAME ${NAME})
    ELSE()
      SET(TPL_IMPORTED_NAME Kokkos::${NAME})
    ENDIF()
  ENDIF()

  IF (NOT TPL_LIBRARY_SUFFIXES)
    SET(TPL_LIBRARY_SUFFIXES lib lib64)
  ENDIF()

  SET(${NAME}_INCLUDE_DIRS)
  IF (TPL_HEADER)
    KOKKOS_FIND_HEADER(${NAME}_INCLUDE_DIRS ${TPL_HEADER} ${NAME} ${ALLOW_PATH_FALLBACK_OPT} PATHS ${TPL_HEADER_PATHS})
  ENDIF()

  FOREACH(HEADER ${TPL_HEADERS})
    KOKKOS_FIND_HEADER(HEADER_FIND_TEMP ${HEADER} ${NAME} ${ALLOW_PATH_FALLBACK_OPT} PATHS ${TPL_HEADER_PATHS})
    IF(HEADER_FIND_TEMP)
      LIST(APPEND ${NAME}_INCLUDE_DIRS ${HEADER_FIND_TEMP})
    ENDIF()
  ENDFOREACH()

  SET(${NAME}_LIBRARY)
  IF(TPL_LIBRARY)
    KOKKOS_FIND_LIBRARY(${NAME}_LIBRARY ${TPL_LIBRARY} ${NAME}
      ${ALLOW_PATH_FALLBACK_OPT}
      PATHS ${TPL_LIBRARY_PATHS}
      SUFFIXES ${TPL_LIBRARY_SUFFIXES})
  ENDIF()

  SET(${NAME}_FOUND_LIBRARIES)
  FOREACH(LIB ${TPL_LIBRARIES})
    KOKKOS_FIND_LIBRARY(${LIB}_LOCATION ${LIB} ${NAME}
      ${ALLOW_PATH_FALLBACK_OPT}
      PATHS ${TPL_LIBRARY_PATHS}
      SUFFIXES ${TPL_LIBRARY_SUFFIXES})
    IF(${LIB}_LOCATION)
      LIST(APPEND ${NAME}_FOUND_LIBRARIES ${${LIB}_LOCATION})
    ELSE()
      SET(${NAME}_FOUND_LIBRARIES ${${LIB}_LOCATION})
      BREAK()
    ENDIF()
  ENDFOREACH()

  INCLUDE(FindPackageHandleStandardArgs)
  #Collect all the variables we need to be valid for
  #find_package to have succeeded
  SET(TPL_VARS_NEEDED)
  IF (TPL_LIBRARY)
    LIST(APPEND TPL_VARS_NEEDED ${NAME}_LIBRARY)
  ENDIF()
  IF(TPL_HEADER)
    LIST(APPEND TPL_VARS_NEEDED ${NAME}_INCLUDE_DIRS)
  ENDIF()
  IF(TPL_LIBRARIES)
    LIST(APPEND TPL_VARS_NEEDED ${NAME}_FOUND_LIBRARIES)
  ENDIF()
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(${TPL_MODULE_NAME} REQUIRED_VARS ${TPL_VARS_NEEDED})

  MARK_AS_ADVANCED(${NAME}_INCLUDE_DIRS ${NAME}_FOUND_LIBRARIES ${NAME}_LIBRARY)

  #this is so much fun on a Cray system
  #/usr/include should never be added as a -isystem include
  #this freaks out the compiler include search order
  IF (KOKKOS_IS_CRAYPE)
    LIST(REMOVE_ITEM ${NAME}_INCLUDE_DIRS "/usr/include")
  ENDIF()

  IF (${TPL_MODULE_NAME}_FOUND)
    SET(IMPORT_TYPE)
    IF (TPL_INTERFACE)
      SET(IMPORT_TYPE "INTERFACE")
    ENDIF()
    KOKKOS_CREATE_IMPORTED_TPL(${TPL_IMPORTED_NAME}
      ${IMPORT_TYPE}
      INCLUDES "${${NAME}_INCLUDE_DIRS}"
      LIBRARY  "${${NAME}_LIBRARY}"
      LINK_LIBRARIES "${${NAME}_FOUND_LIBRARIES}")
  ENDIF()
ENDMACRO(kokkos_find_imported)

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
FUNCTION(kokkos_link_tpl TARGET)
  CMAKE_PARSE_ARGUMENTS(TPL
   "PUBLIC;PRIVATE;INTERFACE"
   "IMPORTED_NAME"
   ""
   ${ARGN})
  #the name of the TPL
  SET(TPL ${TPL_UNPARSED_ARGUMENTS})
  IF (KOKKOS_HAS_TRILINOS)
    #Do nothing, they will have already been linked
  ELSE()
    IF (NOT TPL_IMPORTED_NAME)
      SET(TPL_IMPORTED_NAME Kokkos::${TPL})
    ENDIF()
    IF (KOKKOS_ENABLE_${TPL})
      IF (TPL_PUBLIC)
        TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${TPL_IMPORTED_NAME})
      ELSEIF (TPL_PRIVATE)
        TARGET_LINK_LIBRARIES(${TARGET} PRIVATE ${TPL_IMPORTED_NAME})
      ELSEIF (TPL_INTERFACE)
        TARGET_LINK_LIBRARIES(${TARGET} INTERFACE ${TPL_IMPORTED_NAME})
      ELSE()
        TARGET_LINK_LIBRARIES(${TARGET} ${TPL_IMPORTED_NAME})
      ENDIF()
    ENDIF()
  ENDIF()
ENDFUNCTION()

FUNCTION(COMPILER_SPECIFIC_OPTIONS_HELPER)
  SET(COMPILERS NVIDIA PGI XL DEFAULT Cray Intel Clang AppleClang IntelClang GNU HIP Fujitsu)
  CMAKE_PARSE_ARGUMENTS(
    PARSE
    "LINK_OPTIONS;COMPILE_OPTIONS;COMPILE_DEFINITIONS;LINK_LIBRARIES"
    "COMPILER_ID"
    "${COMPILERS}"
    ${ARGN})
  IF(PARSE_UNPARSED_ARGUMENTS)
    MESSAGE(SEND_ERROR "'${PARSE_UNPARSED_ARGUMENTS}' argument(s) not recognized when providing compiler specific options")
  ENDIF()

  IF(PARSE_COMPILER_ID)
    SET(COMPILER ${${PARSE_COMPILER_ID}})
  ELSE()
    SET(COMPILER ${KOKKOS_CXX_COMPILER_ID})
  ENDIF()

  SET(COMPILER_SPECIFIC_FLAGS_TMP)
  FOREACH(COMP ${COMPILERS})
    IF (COMPILER STREQUAL "${COMP}")
      IF (PARSE_${COMPILER})
        IF (NOT "${PARSE_${COMPILER}}" STREQUAL "NO-VALUE-SPECIFIED")
           SET(COMPILER_SPECIFIC_FLAGS_TMP ${PARSE_${COMPILER}})
        ENDIF()
      ELSEIF(PARSE_DEFAULT)
        SET(COMPILER_SPECIFIC_FLAGS_TMP ${PARSE_DEFAULT})
      ENDIF()
    ENDIF()
  ENDFOREACH()

  IF (PARSE_COMPILE_OPTIONS)
    # The funky logic here is for future handling of argument deduplication
    # If we naively pass multiple -Xcompiler flags to target_compile_options
    # -Xcompiler will get deduplicated and break the build
    IF ("-Xcompiler" IN_LIST COMPILER_SPECIFIC_FLAGS_TMP)
      LIST(REMOVE_ITEM COMPILER_SPECIFIC_FLAGS_TMP "-Xcompiler")
      GLOBAL_APPEND(KOKKOS_XCOMPILER_OPTIONS ${COMPILER_SPECIFIC_FLAGS_TMP})
    ELSE()
      GLOBAL_APPEND(KOKKOS_COMPILE_OPTIONS ${COMPILER_SPECIFIC_FLAGS_TMP})
    ENDIF()
  ENDIF()

  IF (PARSE_LINK_OPTIONS)
    GLOBAL_APPEND(KOKKOS_LINK_OPTIONS ${COMPILER_SPECIFIC_FLAGS_TMP})
  ENDIF()

  IF (PARSE_COMPILE_DEFINITIONS)
    GLOBAL_APPEND(KOKKOS_COMPILE_DEFINITIONS ${COMPILER_SPECIFIC_FLAGS_TMP})
  ENDIF()

  IF (PARSE_LINK_LIBRARIES)
    GLOBAL_APPEND(KOKKOS_LINK_LIBRARIES ${COMPILER_SPECIFIC_FLAGS_TMP})
  ENDIF()
ENDFUNCTION(COMPILER_SPECIFIC_OPTIONS_HELPER)

FUNCTION(COMPILER_SPECIFIC_FLAGS)
  COMPILER_SPECIFIC_OPTIONS_HELPER(${ARGN} COMPILE_OPTIONS LINK_OPTIONS)
ENDFUNCTION(COMPILER_SPECIFIC_FLAGS)

FUNCTION(COMPILER_SPECIFIC_OPTIONS)
  COMPILER_SPECIFIC_OPTIONS_HELPER(${ARGN} COMPILE_OPTIONS)
ENDFUNCTION(COMPILER_SPECIFIC_OPTIONS)

FUNCTION(COMPILER_SPECIFIC_LINK_OPTIONS)
  COMPILER_SPECIFIC_OPTIONS_HELPER(${ARGN} LINK_OPTIONS)
ENDFUNCTION(COMPILER_SPECIFIC_LINK_OPTIONS)

FUNCTION(COMPILER_SPECIFIC_DEFS)
  COMPILER_SPECIFIC_OPTIONS_HELPER(${ARGN} COMPILE_DEFINITIONS)
ENDFUNCTION(COMPILER_SPECIFIC_DEFS)

FUNCTION(COMPILER_SPECIFIC_LIBS)
  COMPILER_SPECIFIC_OPTIONS_HELPER(${ARGN} LINK_LIBRARIES)
ENDFUNCTION(COMPILER_SPECIFIC_LIBS)
# Given a list of the form
#  key1;value1;key2;value2,...
# Create a list of all keys in a variable named ${KEY_LIST_NAME}
# and set the value for each key in a variable ${VAR_PREFIX}key1,...
# kokkos_key_value_map(ARCH ALL_ARCHES key1;value1;key2;value2)
# would produce a list variable ALL_ARCHES=key1;key2
# and individual variables ARCHkey1=value1 and ARCHkey2=value2
MACRO(KOKKOS_KEY_VALUE_MAP VAR_PREFIX KEY_LIST_NAME)
  SET(PARSE_KEY ON)
  SET(${KEY_LIST_NAME})
  FOREACH(ENTRY ${ARGN})
    IF(PARSE_KEY)
      SET(CURRENT_KEY ${ENTRY})
      SET(PARSE_KEY OFF)
      LIST(APPEND ${KEY_LIST_NAME} ${CURRENT_KEY})
    ELSE()
      SET(${VAR_PREFIX}${CURRENT_KEY} ${ENTRY})
      SET(PARSE_KEY ON)
    ENDIF()
  ENDFOREACH()
ENDMACRO()

FUNCTION(KOKKOS_CHECK_DEPRECATED_OPTIONS)
  KOKKOS_KEY_VALUE_MAP(DEPRECATED_MSG_ DEPRECATED_LIST ${ARGN})
  FOREACH(OPTION_SUFFIX ${DEPRECATED_LIST})
    SET(OPTION_NAME Kokkos_${OPTION_SUFFIX})
    SET(OPTION_MESSAGE ${DEPRECATED_MSG_${OPTION_SUFFIX}})
    IF(DEFINED ${OPTION_NAME}) # This variable has been given by the user as on or off
      MESSAGE(SEND_ERROR "Removed option ${OPTION_NAME} has been given with value ${${OPTION_NAME}}. ${OPT_MESSAGE}")
    ENDIF()
  ENDFOREACH()
ENDFUNCTION()

# this function checks whether the current CXX compiler supports building CUDA
FUNCTION(kokkos_cxx_compiler_cuda_test _VAR)
    # don't run this test every time
    IF(DEFINED ${_VAR})
        RETURN()
    ENDIF()

    FILE(WRITE ${PROJECT_BINARY_DIR}/compile_tests/compiles_cuda.cpp
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
")

    TRY_COMPILE(_RET
        ${PROJECT_BINARY_DIR}/compile_tests
        SOURCES ${PROJECT_BINARY_DIR}/compile_tests/compiles_cuda.cpp)

    SET(${_VAR} ${_RET} CACHE STRING "CXX compiler supports building CUDA")
ENDFUNCTION()

# this function is provided to easily select which files use nvcc_wrapper:
#
#       GLOBAL      --> all files
#       TARGET      --> all files in a target
#       SOURCE      --> specific source files
#       DIRECTORY   --> all files in directory
#       PROJECT     --> all files/targets in a project/subproject
#
FUNCTION(kokkos_compilation)
    # check whether the compiler already supports building CUDA
    KOKKOS_CXX_COMPILER_CUDA_TEST(Kokkos_CXX_COMPILER_COMPILES_CUDA)
    # if CUDA compile test has already been performed, just return
    IF(Kokkos_CXX_COMPILER_COMPILES_CUDA)
        RETURN()
    ENDIF()

    CMAKE_PARSE_ARGUMENTS(COMP "GLOBAL;PROJECT" "" "DIRECTORY;TARGET;SOURCE" ${ARGN})

    # find kokkos_launch_compiler
    FIND_PROGRAM(Kokkos_COMPILE_LAUNCHER
        NAMES           kokkos_launch_compiler
        HINTS           ${PROJECT_SOURCE_DIR}
        PATHS           ${PROJECT_SOURCE_DIR}
        PATH_SUFFIXES   bin)

    IF(NOT Kokkos_COMPILE_LAUNCHER)
        MESSAGE(FATAL_ERROR "Kokkos could not find 'kokkos_launch_compiler'. Please set '-DKokkos_COMPILE_LAUNCHER=/path/to/launcher'")
    ENDIF()

    IF(COMP_GLOBAL)
        # if global, don't bother setting others
        SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${Kokkos_COMPILE_LAUNCHER} ${CMAKE_CXX_COMPILER}")
        SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_LINK "${Kokkos_COMPILE_LAUNCHER} ${CMAKE_CXX_COMPILER}")
    ELSE()
        FOREACH(_TYPE PROJECT DIRECTORY TARGET SOURCE)
            # make project/subproject scoping easy, e.g. KokkosCompilation(PROJECT) after project(...)
            IF("${_TYPE}" STREQUAL "PROJECT" AND COMP_${_TYPE})
                LIST(APPEND COMP_DIRECTORY ${PROJECT_SOURCE_DIR})
                UNSET(COMP_${_TYPE})
            ENDIF()
            # set the properties if defined
            IF(COMP_${_TYPE})
                # MESSAGE(STATUS "Using nvcc_wrapper :: ${_TYPE} :: ${COMP_${_TYPE}}")
                SET_PROPERTY(${_TYPE} ${COMP_${_TYPE}} PROPERTY RULE_LAUNCH_COMPILE "${Kokkos_COMPILE_LAUNCHER} ${CMAKE_CXX_COMPILER}")
                SET_PROPERTY(${_TYPE} ${COMP_${_TYPE}} PROPERTY RULE_LAUNCH_LINK "${Kokkos_COMPILE_LAUNCHER} ${CMAKE_CXX_COMPILER}")
            ENDIF()
        ENDFOREACH()
    ENDIF()
ENDFUNCTION()
## KOKKOS_CONFIG_HEADER - parse the data list which is a list of backend names
##                        and create output config header file...used for
##                        creating dynamic include files based on enabled backends
##
##                        SRC_FILE is input file
##                        TARGET_FILE output file
##                        HEADER_GUARD TEXT used with include header guard
##                        HEADER_PREFIX prefix used with include (i.e. fwd, decl, setup)
##                        DATA_LIST list of backends to include in generated file
FUNCTION(KOKKOS_CONFIG_HEADER SRC_FILE TARGET_FILE HEADER_GUARD HEADER_PREFIX DATA_LIST)
   SET(HEADER_GUARD_TAG "${HEADER_GUARD}_HPP_")
   CONFIGURE_FILE(cmake/${SRC_FILE} ${PROJECT_BINARY_DIR}/temp/${TARGET_FILE}.work COPYONLY)
   FOREACH( BACKEND_NAME ${DATA_LIST} )
   SET(INCLUDE_NEXT_FILE "#include <${HEADER_PREFIX}_${BACKEND_NAME}.hpp>
\@INCLUDE_NEXT_FILE\@")
   CONFIGURE_FILE(${PROJECT_BINARY_DIR}/temp/${TARGET_FILE}.work ${PROJECT_BINARY_DIR}/temp/${TARGET_FILE}.work @ONLY)
   ENDFOREACH()
   SET(INCLUDE_NEXT_FILE "" )
   CONFIGURE_FILE(${PROJECT_BINARY_DIR}/temp/${TARGET_FILE}.work ${TARGET_FILE} @ONLY)
ENDFUNCTION()
