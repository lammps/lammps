#FIXME_OPENMPTARGET - compiling in debug mode causes ICE.
KOKKOS_ADD_BENCHMARK_DIRECTORIES(atomic)
KOKKOS_ADD_BENCHMARK_DIRECTORIES(gather)
KOKKOS_ADD_BENCHMARK_DIRECTORIES(gups)
KOKKOS_ADD_BENCHMARK_DIRECTORIES(launch_latency)
KOKKOS_ADD_BENCHMARK_DIRECTORIES(stream)
KOKKOS_ADD_BENCHMARK_DIRECTORIES(view_copy_constructor)
#FIXME_OPENMPTARGET - These two benchmarks cause ICE. Commenting them for now but a deeper analysis on the cause and a possible fix will follow.
IF(NOT Kokkos_ENABLE_OPENMPTARGET)
    KOKKOS_ADD_BENCHMARK_DIRECTORIES(policy_performance)
    KOKKOS_ADD_BENCHMARK_DIRECTORIES(bytes_and_flops)
ENDIF()
