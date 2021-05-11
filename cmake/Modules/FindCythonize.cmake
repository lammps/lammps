# Find the Cythonize tool.
#
# This code sets the following variables:
#
#  Cythonize_EXECUTABLE
#
# adapted from https://github.com/cmarshall108/cython-cmake-example/blob/master/cmake/FindCython.cmake
#=============================================================================

if(CMAKE_VERSION VERSION_LESS 3.12)
    find_package(PythonInterp 3.6 QUIET) # Deprecated since version 3.12
    if(PYTHONINTERP_FOUND)
        set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
    endif()
else()
    find_package(Python3 3.6 COMPONENTS Interpreter QUIET)
endif()

# Use the Cython executable that lives next to the Python executable
# if it is a local installation.
if(Python3_EXECUTABLE)
  get_filename_component(_python_path ${Python3_EXECUTABLE} PATH)
  find_program(Cythonize_EXECUTABLE
    NAMES cythonize3 cythonize cythonize.bat
    HINTS ${_python_path})
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Cythonize REQUIRED_VARS Cythonize_EXECUTABLE)
mark_as_advanced(Cythonize_EXECUTABLE)
