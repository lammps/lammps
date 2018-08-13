if (KOKKOS_ENABLE_CUDA)
  add_subdirectory(Cuda)
endif()
if (KOKKOS_ENABLE_OPENMP)
  add_subdirectory(OpenMP)
endif()
if (KOKKOS_ENABLE_ROCM)
  add_subdirectory(ROCm)
endif()
if (KOKKOS_ENABLE_SERIAL)
  add_subdirectory(Serial)
endif()
if (KOKKOS_ENABLE_THREADS)
  add_subdirectory(Threads)
endif()

set(ETI_SOURCES "${ETI_SOURCES}" PARENT_SCOPE)

install(FILES
common/Kokkos_ViewFillCopyETIAvail_Macros.hpp
common/Kokkos_ViewFillCopyETIDecl_Macros.hpp
DESTINATION include/eti/common)
