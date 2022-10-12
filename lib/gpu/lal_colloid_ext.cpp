/***************************************************************************
                                 colloid_ext.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to colloid acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_colloid.h"

using namespace std;
using namespace LAMMPS_AL;

static Colloid<PRECISION,ACC_PRECISION> COLLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int colloid_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                     double **host_lj2, double **host_lj3, double **host_lj4,
                     double **offset, double *special_lj,
                     double **host_a12, double **host_a1, double **host_a2,
                     double **host_d1, double **host_d2, double **host_sigma3,
                     double **host_sigma6, int **host_form, const int inum,
                     const int nall, const int max_nbors,  const int maxspecial,
                     const double cell_size, int &gpu_mode, FILE *screen) {
  COLLMF.clear();
  gpu_mode=COLLMF.device->gpu_mode();
  double gpu_split=COLLMF.device->particle_split();
  int first_gpu=COLLMF.device->first_device();
  int last_gpu=COLLMF.device->last_device();
  int world_me=COLLMF.device->world_me();
  int gpu_rank=COLLMF.device->gpu_rank();
  int procs_per_gpu=COLLMF.device->procs_per_gpu();

  COLLMF.device->init_message(screen,"colloid",first_gpu,last_gpu);

  bool message=false;
  if (COLLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=COLLMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3,
                        host_lj4, offset, special_lj, host_a12, host_a1,
                        host_a2, host_d1, host_d2, host_sigma3,
                        host_sigma6, host_form, inum, nall, max_nbors,
                        maxspecial, cell_size, gpu_split, screen);

  COLLMF.device->world_barrier();
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
      init_ok=COLLMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                          offset, special_lj, host_a12, host_a1, host_a2,
                          host_d1, host_d2, host_sigma3, host_sigma6, host_form,
                          inum, nall, max_nbors, maxspecial,
                          cell_size, gpu_split, screen);

    COLLMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    COLLMF.estimate_gpu_overhead();
  return init_ok;
}

void colloid_gpu_clear() {
  COLLMF.clear();
}

int ** colloid_gpu_compute_n(const int ago, const int inum_full,
                             const int nall, double **host_x, int *host_type,
                             double *sublo, double *subhi, tagint *tag, int **nspecial,
                             tagint **special, const bool eflag, const bool vflag,
                             const bool eatom, const bool vatom, int &host_start,
                             int **ilist, int **jnum, const double cpu_time,
                             bool &success) {
  return COLLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, ilist, jnum, cpu_time, success);
}

void colloid_gpu_compute(const int ago, const int inum_full, const int nall,
                         double **host_x, int *host_type, int *ilist, int *numj,
                         int **firstneigh, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         const double cpu_time, bool &success) {
  COLLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                 firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double colloid_gpu_bytes() {
  return COLLMF.host_memory_usage();
}


