#undef PACKAGE
#define PACKAGE "REPLICA"
#include "REPLICA/verlet_split.h"
#undef PACKAGE
#define PACKAGE "USER-OMP"
#include "USER-OMP/respa_omp.h"
#undef PACKAGE
#define PACKAGE "KOKKOS"
#include "KOKKOS/verlet_kokkos.h"
#undef PACKAGE
#define PACKAGE "USER-INTEL"
#include "USER-INTEL/verlet_lrt_intel.h"
