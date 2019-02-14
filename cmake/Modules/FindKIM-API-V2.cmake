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
# Copyright (c) 2013--2018, Regents of the University of Minnesota.
# All rights reserved.
#
# Contributors:
#    Richard Berger
#    Christoph Junghans
#    Ryan S. Elliott
#

# - Find KIM-API-V2
#
# sets standard pkg_check_modules variables plus:
#
# KIM-API-V2-CMAKE_C_COMPILER
# KIM-API-V2-CMAKE_CXX_COMPILER
# KIM-API-V2-CMAKE_Fortran_COMPILER
#
find_package(PkgConfig REQUIRED)
include(FindPackageHandleStandardArgs)

pkg_check_modules(KIM-API-V2 REQUIRED libkim-api-v2>=2.0)

pkg_get_variable(KIM-API-V2-CMAKE_C_COMPILER libkim-api-v2 CMAKE_C_COMPILER)
pkg_get_variable(KIM-API-V2-CMAKE_CXX_COMPILER libkim-api-v2 CMAKE_CXX_COMPILER)
pkg_get_variable(KIM-API-V2_CMAKE_Fortran_COMPILER libkim-api-v2 CMAKE_Fortran_COMPILER)

# handle the QUIETLY and REQUIRED arguments and set KIM-API-V2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(KIM-API-V2 REQUIRED_VARS KIM-API-V2_LIBRARIES)
