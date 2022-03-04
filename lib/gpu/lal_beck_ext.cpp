/***************************************************************************
                                 beck_ext.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to beck acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_beck.h"

using namespace std;
using namespace LAMMPS_AL;

static Beck<PRECISION,ACC_PRECISION> BLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int beck_gpu_init(const int ntypes, double **cutsq, double **aa,
                  double **alpha, double **beta, double **AA, double **BB,
                  double *special_lj, const int inum, const int nall,
                  const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen) {
  BLMF.clear();
  gpu_mode=BLMF.device->gpu_mode();
  double gpu_split=BLMF.device->particle_split();
  int first_gpu=BLMF.device->first_device();
  int last_gpu=BLMF.device->last_device();
  int world_me=BLMF.device->world_me();
  int gpu_rank=BLMF.device->gpu_rank();
  int procs_per_gpu=BLMF.device->procs_per_gpu();

  BLMF.device->init_message(screen,"beck",first_gpu,last_gpu);

  bool message=false;
  if (BLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=BLMF.init(ntypes, cutsq, aa, alpha, beta,
                      AA, BB, special_lj, inum, nall, max_nbors,
                      maxspecial, cell_size, gpu_split, screen);

  BLMF.device->world_barrier();
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
      init_ok=BLMF.init(ntypes, cutsq, aa, alpha, beta, AA, BB,
                        special_lj, inum, nall, max_nbors, maxspecial,
                        cell_size, gpu_split, screen);

    BLMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    BLMF.estimate_gpu_overhead();
  return init_ok;
}

void beck_gpu_clear() {
  BLMF.clear();
}

int ** beck_gpu_compute_n(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success) {
  return BLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                      subhi, tag, nspecial, special, eflag, vflag, eatom,
                      vatom, host_start, ilist, jnum, cpu_time, success);
}

void beck_gpu_compute(const int ago, const int inum_full, const int nall,
                       double **host_x, int *host_type, int *ilist, int *numj,
                       int **firstneigh, const bool eflag, const bool vflag,
                       const bool eatom, const bool vatom, int &host_start,
                       const double cpu_time, bool &success) {
  BLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
               firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double beck_gpu_bytes() {
  return BLMF.host_memory_usage();
}


