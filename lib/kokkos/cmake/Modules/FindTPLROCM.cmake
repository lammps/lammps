include(FindPackageHandleStandardArgs)

FIND_LIBRARY(AMD_HIP_LIBRARY amdhip64 PATHS ENV ROCM_PATH PATH_SUFFIXES lib)
FIND_LIBRARY(HSA_RUNTIME_LIBRARY hsa-runtime64 PATHS ENV ROCM_PATH PATH_SUFFIXES lib)

# FIXME_HIP Starting with ROCm 5.5 it is not necessary to link againt clang_rt.
# We keep the code as is for now because it is hard to find the version of ROCM
# found.
# clang_rt.builtins is necessary to use half precision. The following code to
# find clang_rt.buitins is based on
# https://github.com/ROCm-Developer-Tools/hipamd/blob/d1e0ee98a0f3d79f7bf43295f82d0053a69ec742/hip-config.cmake.in#L241
# NOTE: Per the above, we still search for the clang-rt library,
# but use the user's specified compiler to find the library to avoid use of
# environment variables / relative paths.
execute_process(
  COMMAND ${CMAKE_CXX_COMPILER} -print-libgcc-file-name --rtlib=compiler-rt
  OUTPUT_VARIABLE CLANG_RT_LIBRARY
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE CLANG_RT_CHECK)

if( NOT "${CLANG_RT_CHECK}" STREQUAL "0" )
  # if the above failed, we delete CLANG_RT_LIBRARY to make the args check
  # below fail
  unset(CLANG_RT_LIBRARY)
endif()


find_package_handle_standard_args(TPLROCM DEFAULT_MSG AMD_HIP_LIBRARY HSA_RUNTIME_LIBRARY CLANG_RT_LIBRARY)

kokkos_create_imported_tpl(ROCM INTERFACE
  LINK_LIBRARIES ${HSA_RUNTIME_LIBRARY} ${AMD_HIP_LIBRARY} ${CLANG_RT_LIBRARY}
  COMPILE_DEFINITIONS __HIP_ROCclr__
)
