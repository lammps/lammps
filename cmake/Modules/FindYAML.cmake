# - Find libyaml
# Find the native Yaml headers and libraries.
#
#  YAML_INCLUDE_DIRS - where to find yaml.h
#  YAML_LIBRARIES    - List of libraries when using libyaml
#  YAML_FOUND        - True if libyaml is found.
#

find_path(YAML_INCLUDE_DIR yaml.h PATH_SUFFIXES yaml)
find_library(YAML_LIBRARY NAMES yaml)

# handle the QUIET and REQUIRED arguments and
# set YAML_FOUND to TRUE if all variables are non-zero
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YAML DEFAULT_MSG YAML_LIBRARY YAML_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(YAML_FOUND)
  set(YAML_LIBRARIES ${YAML_LIBRARY})
  set(YAML_INCLUDE_DIRS ${YAML_INCLUDE_DIR})

  if(NOT TARGET Yaml::Yaml)
    add_library(Yaml::Yaml UNKNOWN IMPORTED)
    set_target_properties(Yaml::Yaml PROPERTIES
      IMPORTED_LOCATION "${YAML_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${YAML_INCLUDE_DIR}")
  endif()
endif()

mark_as_advanced(YAML_INCLUDE_DIR YAML_LIBRARY)
