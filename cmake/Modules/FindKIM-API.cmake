#
# CDDL HEADER START
#
# The contents of this file are subject to the terms of the Common Development
# and Distribution License Version 1.0 (the "License").
#
# You can obtain a copy of the license at
# http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
# specific language governing permissions and limitations under the License.
#
# When distributing Covered Code, include this CDDL HEADER in each file and
# include the License file in a prominent location with the name LICENSE.CDDL.
# If applicable, add the following below this CDDL HEADER, with the fields
# enclosed by brackets "[]" replaced with your own identifying information:
#
# Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
#
# CDDL HEADER END
#

#
# Copyright (c) 2013--2019, Regents of the University of Minnesota.
# All rights reserved.
#
# Contributors:
#    Richard Berger
#    Christoph Junghans
#    Ryan S. Elliott
#

# - Find KIM-API
#
# sets standard pkg_check_modules variables plus:
#
# KIM-API-CMAKE_C_COMPILER
# KIM-API-CMAKE_CXX_COMPILER
# KIM-API-CMAKE_Fortran_COMPILER
#

function(_KIMAPI_GET_VERSION _OUT_ver _version_hdr)
  if(NOT EXISTS ${_version_hdr})
    message(FATAL_ERROR "Header file ${_version_hdr} not found (check value of KIM-API_INCLUDE_DIR)")
  endif()
  foreach(_var KIM_VERSION_MAJOR KIM_VERSION_MINOR KIM_VERSION_PATCH)  
    file(STRINGS ${_version_hdr} _contents REGEX "#define ${_var}[ \t]+")
    if(_contents)
      string(REGEX REPLACE ".*#define ${_var}[ \t]+([0-9]+).*" "\\1" _${_var} "${_contents}")
      if(${_${_var}} STREQUAL "")
        message(FATAL_ERROR "Version parsing failed for ${_var} in ${_version_hdr}, got empty return!")
      elseif(NOT ${_${_var}} MATCHES "^[0-9]+$")
        message(FATAL_ERROR "Version parsing failed for ${_var} in ${_version_hdr}, excepted a number but got ${_${_var}}!")
      endif()
    else()
      message(FATAL_ERROR "No ${_var} line found in include file ${_version_hdr}")
    endif()
  endforeach()
  set(${_OUT_ver} ${_KIM_VERSION_MAJOR}.${_KIM_VERSION_MINOR}.${_KIM_VERSION_PATCH} PARENT_SCOPE)
endfunction()

if(KIM-API_FIND_QUIETLY)
  set(REQ_OR_QUI "QUIET")
else()
  set(REQ_OR_QUI "REQUIRED")
endif()

find_package(PkgConfig ${REQ_OR_QUI})
include(FindPackageHandleStandardArgs)

pkg_check_modules(KIM-API ${REQ_OR_QUI} libkim-api>=2.0)

if(KIM-API_FOUND)
  pkg_get_variable(KIM-API-CMAKE_C_COMPILER libkim-api CMAKE_C_COMPILER)
  pkg_get_variable(KIM-API-CMAKE_CXX_COMPILER libkim-api CMAKE_CXX_COMPILER)
  pkg_get_variable(KIM-API_CMAKE_Fortran_COMPILER libkim-api CMAKE_Fortran_COMPILER)
endif()

if(KIM-API_INCLUDEDIR)
  _KIMAPI_GET_VERSION(KIM-API_VERSION ${KIM-API_INCLUDEDIR}/KIM_Version.h)
else()
  message(FATAL_ERROR "KIM-API_INCLUDEDIR not defined")
endif()

# handle the QUIETLY and REQUIRED arguments and set KIM-API_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(KIM-API REQUIRED_VARS KIM-API_LIBRARIES VERSION_VAR KIM-API_VERSION)
