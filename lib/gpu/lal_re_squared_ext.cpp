/***************************************************************************
                              re_squared_ext.cpp
                             -------------------
                               W. Michael Brown

  LAMMPS Wrappers for RE-Squared Acceleration

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_re_squared.h"

using namespace std;
using namespace LAMMPS_AL;

static RESquared<PRECISION,ACC_PRECISION> REMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int re_gpu_init(const int ntypes, double **shape, double **well, double **cutsq,
                double **sigma, double **epsilon,
                int **form, double **host_lj1, double **host_lj2,
                double **host_lj3, double **host_lj4, double **offset,
                double *special_lj, const int inum, const int nall,
                const int max_nbors, const int maxspecial,
                const double cell_size, int &gpu_mode, FILE *screen) {
  REMF.clear();
  gpu_mode=REMF.device->gpu_mode();
  double gpu_split=REMF.device->particle_split();
  int first_gpu=REMF.device->first_device();
  int last_gpu=REMF.device->last_device();
  int world_me=REMF.device->world_me();
  int gpu_rank=REMF.device->gpu_rank();
  int procs_per_gpu=REMF.device->procs_per_gpu();

  REMF.device->init_message(screen,"resquared",first_gpu,last_gpu);

  bool message=false;
  if (REMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=REMF.init(ntypes, shape, well, cutsq, sigma, epsilon,
                      form, host_lj1, host_lj2, host_lj3, host_lj4, offset,
                      special_lj, inum, nall, max_nbors, maxspecial, cell_size,
                      gpu_split, screen);

  REMF.device->world_barrier();
  if (message)
    fprintf(screen,"Done.\n");

  for (int i=0; i<procs_per_gpu; i++) {
    if (message) {
      if (last_gpu-first_gpu==0)
        fprintf(screen,"Initializing Device %d on core %d...",first_gpu,i);
      else
        fprintf(screen,"Initializing Devices %d-%d on core %d...",first_gpu,
                last_gpu,i);
      fflush(screen);
    }
    if (gpu_rank==i && world_me!=0)
      init_ok=REMF.init(ntypes, shape, well, cutsq,  sigma, epsilon,
                        form, host_lj1, host_lj2, host_lj3,
                        host_lj4, offset, special_lj,  inum, nall,
                        max_nbors, maxspecial, cell_size, gpu_split, screen);

    REMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    REMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Clear memory on host and device
// ---------------------------------------------------------------------------
void re_gpu_clear() {
  REMF.clear();
}

  int** compute(const int ago, const int inum_full, const int nall,
                double **host_x, int *host_type, double *sublo,
                double *subhi, tagint *tag, int **nspecial,
                tagint **special, const bool eflag, const bool vflag,
                const bool eatom, const bool vatom, int &host_start,
                int **ilist, int **numj, const double cpu_time, bool &success,
                double **host_quat);

int** re_gpu_compute_n(const int ago, const int inum_full, const int nall,
                       double **host_x, int *host_type, double *sublo,
                       double *subhi, tagint *tag, int **nspecial, tagint **special,
                       const bool eflag, const bool vflag, const bool eatom,
                       const bool vatom, int &host_start, int **ilist,
                       int **jnum, const double cpu_time, bool &success,
                       double **host_quat) {
  return REMF.compute(ago, inum_full, nall, host_x, host_type, sublo, subhi,
                      tag, nspecial, special, eflag, vflag, eatom, vatom,
                      host_start, ilist, jnum, cpu_time, success, host_quat);
}

int * re_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success, double **host_quat) {
  return REMF.compute(ago, inum_full, nall, host_x, host_type, ilist,
                      numj, firstneigh, eflag, vflag, eatom, vatom, host_start,
                      cpu_time, success, host_quat);
}

// ---------------------------------------------------------------------------
// Return memory usage
// ---------------------------------------------------------------------------
double re_gpu_bytes() {
  return REMF.host_memory_usage();
}

