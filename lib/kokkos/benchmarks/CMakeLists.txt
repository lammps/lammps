#FIXME_OPENMPTARGET - compiling in debug mode causes ICE.
kokkos_add_benchmark_directories(atomic)
kokkos_add_benchmark_directories(gather)
kokkos_add_benchmark_directories(gups)
kokkos_add_benchmark_directories(launch_latency)
kokkos_add_benchmark_directories(stream)
kokkos_add_benchmark_directories(view_copy_constructor)
#FIXME_OPENMPTARGET - These two benchmarks cause ICE. Commenting them for now but a deeper analysis on the cause and a possible fix will follow.
if(NOT Kokkos_ENABLE_OPENMPTARGET)
  kokkos_add_benchmark_directories(policy_performance)
  kokkos_add_benchmark_directories(bytes_and_flops)
endif()
