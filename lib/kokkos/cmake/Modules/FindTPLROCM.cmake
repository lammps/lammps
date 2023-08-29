include(FindPackageHandleStandardArgs)

FIND_LIBRARY(AMD_HIP_LIBRARY amdhip64 PATHS ENV ROCM_PATH PATH_SUFFIXES lib)
FIND_LIBRARY(HSA_RUNTIME_LIBRARY hsa-runtime64 PATHS ENV ROCM_PATH PATH_SUFFIXES lib)

# clang_rt.builtins is necessary to use half precision. The following code to
# find clang_rt.buitins is based on
# https://github.com/ROCm-Developer-Tools/HIP/blob/develop/hip-lang-config.cmake.in#L99-L111
file(GLOB_RECURSE CLANG_RT_DIR "$ENV{ROCM_PATH}/llvm/lib/clang/*/lib/*/*clang_rt.builtins*")
FIND_LIBRARY(CLANG_RT_LIBRARY
  NAMES
  clang_rt.builtins
  clang_rt.builtins-x86_64
  PATHS "${CLANG_RT_DIR}/..")

find_package_handle_standard_args(TPLROCM DEFAULT_MSG AMD_HIP_LIBRARY HSA_RUNTIME_LIBRARY CLANG_RT_LIBRARY)

kokkos_create_imported_tpl(ROCM INTERFACE
  LINK_LIBRARIES ${HSA_RUNTIME_LIBRARY} ${AMD_HIP_LIBRARY} ${CLANG_RT_LIBRARY}
  COMPILE_DEFINITIONS __HIP_ROCclr__
)
