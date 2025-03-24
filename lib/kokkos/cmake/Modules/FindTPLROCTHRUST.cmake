# ROCm 5.6 and earlier set AMDGPU_TARGETS and GPU_TARGETS to all the supported
# architectures. Therefore, we end up compiling Kokkos for all the supported
# architecture. Starting with ROCm 5.7 AMDGPU_TARGETS and GPU_TARGETS are empty.
# It is the user's job to set the variables. Since we are injecting the
# architecture flag ourselves, we can let the variables empty. To replicate the
# behavior of ROCm 5.7 and later for earlier version of ROCm we set
# AMDGPU_TARGETS and GPU_TARGETS to empty and set the values in the cache. If
# the values are not cached, FIND_PACKAGE(rocthrust) will overwrite them.
set(AMDGPU_TARGETS "" CACHE STRING "AMD GPU targets to compile for")
set(GPU_TARGETS "" CACHE STRING "GPU targets to compile for")
find_package(rocthrust REQUIRED)
kokkos_create_imported_tpl(ROCTHRUST INTERFACE LINK_LIBRARIES roc::rocthrust)

# Export ROCTHRUST as a Kokkos dependency
kokkos_export_cmake_tpl(rocthrust)
