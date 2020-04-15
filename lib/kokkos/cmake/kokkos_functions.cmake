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
    GET_TARGET_PROPERTY(LIB_TYPE ${NAME} TYPE)
    IF (${LIB_TYPE} STREQUAL "INTERFACE_LIBRARY")
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
      KOKKOS_APPEND_CONFIG_LINE("ADD_LIBRARY(${NAME} UNKNOWN IMPORTED)")
      KOKKOS_APPEND_CONFIG_LINE("SET_TARGET_PROPERTIES(${NAME} PROPERTIES")
      
      GET_TARGET_PROPERTY(TPL_LIBRARY ${NAME} IMPORTED_LOCATION)
      IF(TPL_LIBRARY)
        KOKKOS_APPEND_CONFIG_LINE("IMPORTED_LOCATION ${TPL_LIBRARY}")
      ENDIF()

      GET_TARGET_PROPERTY(TPL_INCLUDES ${NAME} INTERFACE_INCLUDE_DIRECTORIES)
      IF(TPL_INCLUDES)
        KOKKOS_APPEND_CONFIG_LINE("INTERFACE_INCLUDE_DIRECTORIES ${TPL_INCLUDES}")
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
        KOKKOS_APPEND_CONFIG_LINE("INTERFACE_LINK_LIBRARIES ${TPL_LINK_LIBRARIES}")
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
# 2. Kokkos_<NAME>_DIR variable
# 3. Locations in the PATHS option
# 4. Default system paths, if allowed.
#
# Default system paths are allowed if none of options (1)-(3) are specified
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

  SET(${HEADER}_FOUND FALSE)
  SET(HAVE_CUSTOM_PATHS FALSE)
  IF(NOT ${HEADER}_FOUND AND DEFINED ${TPL_NAME}_ROOT)
    #ONLY look in the root directory
    FIND_PATH(${VAR_NAME} ${HEADER} PATHS ${${TPL_NAME}_ROOT}/include NO_DEFAULT_PATH)
    SET(HAVE_CUSTOM_PATHS TRUE)
  ENDIF()

  IF(NOT ${HEADER}_FOUND AND DEFINED KOKKOS_${TPL_NAME}_DIR)
    #ONLY look in the root directory
    FIND_PATH(${VAR_NAME} ${HEADER} PATHS ${KOKKOS_${TPL_NAME}_DIR}/include NO_DEFAULT_PATH)
    SET(HAVE_CUSTOM_PATHS TRUE)
  ENDIF()

  IF (NOT ${HEADER}_FOUND AND TPL_PATHS)
    #we got custom paths
    #ONLY look in these paths and nowhere else
    FIND_PATH(${VAR_NAME} ${HEADER} PATHS ${TPL_PATHS} NO_DEFAULT_PATH)
    SET(HAVE_CUSTOM_PATHS TRUE)
  ENDIF()

  IF (NOT HAVE_CUSTOM_PATHS OR TPL_ALLOW_SYSTEM_PATH_FALLBACK)
    #Now go ahead and look in system paths
    IF (NOT ${HEADER}_FOUND)
      FIND_PATH(${VAR_NAME} ${HEADER})
    ENDIF()
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
# 2. Kokkos_<NAME>_DIR variable
# 3. Locations in the PATHS option
# 4. Default system paths, if allowed.
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
MACRO(kokkos_find_library VAR_NAME LIB TPL_NAME)
  CMAKE_PARSE_ARGUMENTS(TPL
   "ALLOW_SYSTEM_PATH_FALLBACK"
   ""
   "PATHS"
   ${ARGN})

  SET(${LIB}_FOUND FALSE)
  SET(HAVE_CUSTOM_PATHS FALSE)
  IF(NOT ${LIB}_FOUND AND DEFINED ${TPL_NAME}_ROOT)
    FIND_LIBRARY(${VAR_NAME} ${LIB} PATHS ${${TPL_NAME}_ROOT}/lib ${${TPL_NAME}_ROOT}/lib64 NO_DEFAULT_PATH)
    SET(HAVE_CUSTOM_PATHS TRUE)
  ENDIF()

  IF(NOT ${LIB}_FOUND AND DEFINED KOKKOS_${TPL_NAME}_DIR)
    #we got root paths, only look in these paths and nowhere else
    FIND_LIBRARY(${VAR_NAME} ${LIB} PATHS ${KOKKOS_${TPL_NAME}_DIR}/lib ${KOKKOS_${TPL_NAME}_DIR}/lib64 NO_DEFAULT_PATH)
    SET(HAVE_CUSTOM_PATHS TRUE)
  ENDIF()

  IF (NOT ${LIB}_FOUND AND TPL_PATHS)
    #we got custom paths, only look in these paths and nowhere else
    FIND_LIBRARY(${VAR_NAME} ${LIB} PATHS ${TPL_PATHS} NO_DEFAULT_PATH)
    SET(HAVE_CUSTOM_PATHS TRUE)
  ENDIF()


  IF (NOT HAVE_CUSTOM_PATHS OR TPL_ALLOW_SYSTEM_PATH_FALLBACK)
    IF (NOT ${LIB}_FOUND)
      #Now go ahead and look in system paths
      FIND_LIBRARY(${VAR_NAME} ${LIB})
    ENDIF()
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
#   ``ALLOW_SYSTEM_PATH_FALLBACK"
#
#     If custom paths are given and the library is not found
#     should we be allowed to search default system paths
#     or error out if not found in given paths.
#
#   ``LIBRARY <name>``
#
#     If specified, this gives the name of the library to look for
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
#   ``LIBRARY_PATHS <path1> <path2> ...``
#
#     If specified, this gives a list of paths to search for the library
#     If not given, <NAME>_ROOT/lib and <NAME>_ROOT/lib64 will be searched.
#
#   ``HEADER_PATHS <path1> <path2> ...``
#
#     If specified, this gives a list of paths to search for the headers
#     If not given, <NAME>_ROOT/include and <NAME>_ROOT/include will be searched.
#
#   ``HEADERS <name1> <name2> ...``
#
#     If specified, this gives a list of headers to find for the package
#
#   ``LIBRARIES <name1> <name2> ...``
#
#     If specified, this gives a list of libraries to find for the package
#
MACRO(kokkos_find_imported NAME)
  CMAKE_PARSE_ARGUMENTS(TPL
   "INTERFACE;ALLOW_SYSTEM_PATH_FALLBACK"
   "HEADER;LIBRARY;IMPORTED_NAME;MODULE_NAME"
   "HEADER_PATHS;LIBRARY_PATHS;HEADERS;LIBRARIES"
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
    KOKKOS_FIND_LIBRARY(${NAME}_LIBRARY ${TPL_LIBRARY} ${NAME} ${ALLOW_PATH_FALLBACK_OPT} PATHS ${TPL_LIBRARY_PATHS})
  ENDIF()

  SET(${NAME}_FOUND_LIBRARIES)
  FOREACH(LIB ${TPL_LIBRARIES})
    KOKKOS_FIND_LIBRARY(${LIB}_LOCATION ${LIB} ${NAME} ${ALLOW_PATH_FALLBACK_OPT} PATHS ${TPL_LIBRARY_PATHS})
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

