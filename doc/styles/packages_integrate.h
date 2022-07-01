#undef PACKAGE
#define PACKAGE "REPLICA"
#include "REPLICA/verlet_split.h"
#undef PACKAGE
#define PACKAGE "KOKKOS"
#include "KOKKOS/verlet_kokkos.h"
#undef PACKAGE
#define PACKAGE "INTEL"
#include "INTEL/verlet_lrt_intel.h"
#undef PACKAGE
#define PACKAGE "OPENMP"
#include "OPENMP/respa_omp.h"
