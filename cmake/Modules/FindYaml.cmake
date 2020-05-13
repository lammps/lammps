# - Find libyaml
# Find the native Yaml headers and libraries.
#
#  Yaml_INCLUDE_DIRS - where to find yaml.h
#  Yaml_LIBRARIES    - List of libraries when using libyaml
#  Yaml_FOUND        - True if libyaml is found.
#

find_path(Yaml_INCLUDE_DIR yaml.h PATH_SUFFIXES yaml)
find_library(Yaml_LIBRARY NAMES yaml)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set YAML_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(Yaml DEFAULT_MSG YAML_LIBRARY YAML_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(Yaml_FOUND)
  set(Yaml_LIBRARIES ${Yaml_LIBRARY})
  set(Yaml_INCLUDE_DIRS ${Yaml_INCLUDE_DIR})

  if(NOT TARGET Yaml::Yaml)
    add_library(Yaml::Yaml UNKNOWN IMPORTED)
    set_target_properties(Yaml::Yaml PROPERTIES
      IMPORTED_LOCATION "${Yaml_LIBRARY}" 
      INTERFACE_INCLUDE_DIRECTORIES "${Yaml_INCLUDE_DIR}")
  endif()
endif()

mark_as_advanced(Yaml_INCLUDE_DIR Yaml_LIBRARY )
