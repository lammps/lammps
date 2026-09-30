find_package(rocthrust REQUIRED PATHS ${ROCM_PATH} $ENV{ROCM_PATH})
kokkos_create_imported_tpl(ROCTHRUST INTERFACE LINK_LIBRARIES roc::rocthrust)

# Export ROCTHRUST as a Kokkos dependency
kokkos_export_cmake_tpl(rocthrust)
