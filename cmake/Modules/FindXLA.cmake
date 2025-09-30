# - Find XLA
# Find the XLA FFI headers and JAX installation.
#
#  XLA_INCLUDE_DIRS - where to find XLA FFI headers
#  XLA_FOUND        - True if XLA/JAX found and usable
#

# Check if JAX is installed in the current Python environment
execute_process(
  COMMAND python -c "import importlib.util; print('YES' if importlib.util.find_spec('jax') else 'NO')"
  OUTPUT_VARIABLE JAX_INSTALLED
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE JAX_CHECK_RESULT)

# Get XLA FFI include directory if JAX is available
if(JAX_CHECK_RESULT EQUAL 0 AND JAX_INSTALLED STREQUAL "YES")
  # This import is needed for the JAX, XLA FFI bindings.
  # More infos at https://docs.jax.dev/en/latest/ffi.html
  # For the following code, having jax cpu version installed is sufficient
  # $ pip install jax
  execute_process(
    COMMAND python -c "from jax import ffi; print(ffi.include_dir())"
    OUTPUT_VARIABLE XLA_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE XLA_FFI_RESULT)
  
  if(XLA_FFI_RESULT EQUAL 0 AND XLA_INCLUDE_DIR)
    set(XLA_FFI_AVAILABLE TRUE)
  else()
    set(XLA_FFI_AVAILABLE FALSE)
  endif()
else()
  set(XLA_FFI_AVAILABLE FALSE)
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set XLA_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(XLA DEFAULT_MSG XLA_FFI_AVAILABLE XLA_INCLUDE_DIR)

if(XLA_FOUND)
  set(XLA_INCLUDE_DIRS ${XLA_INCLUDE_DIR})
  
  if(NOT TARGET XLA::XLA)
    add_library(XLA::XLA INTERFACE IMPORTED)
    set_target_properties(XLA::XLA PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${XLA_INCLUDE_DIR}"
      INTERFACE_COMPILE_DEFINITIONS "WITH_JAX")
  endif()
endif()

mark_as_advanced(XLA_INCLUDE_DIR XLA_FFI_AVAILABLE) 