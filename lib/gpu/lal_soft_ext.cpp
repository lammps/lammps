/***************************************************************************
                                 soft_ext.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to soft acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_soft.h"

using namespace std;
using namespace LAMMPS_AL;

static Soft<PRECISION,ACC_PRECISION> SLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int soft_gpu_init(const int ntypes, double **cutsq, double **host_prefactor,
                  double **host_cut, double *special_lj,
                  const int inum, const int nall, const int max_nbors,
                  const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen) {
  SLMF.clear();
  gpu_mode=SLMF.device->gpu_mode();
  double gpu_split=SLMF.device->particle_split();
  int first_gpu=SLMF.device->first_device();
  int last_gpu=SLMF.device->last_device();
  int world_me=SLMF.device->world_me();
  int gpu_rank=SLMF.device->gpu_rank();
  int procs_per_gpu=SLMF.device->procs_per_gpu();

  SLMF.device->init_message(screen,"soft",first_gpu,last_gpu);

  bool message=false;
  if (SLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=SLMF.init(ntypes, cutsq, host_prefactor, host_cut,
                      special_lj, inum, nall, 300,
                      maxspecial, cell_size, gpu_split, screen);

  SLMF.device->world_barrier();
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
      init_ok=SLMF.init(ntypes, cutsq, host_prefactor, host_cut,
                        special_lj, inum, nall, 300, maxspecial,
                        cell_size, gpu_split, screen);

    SLMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    SLMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated constants to device
// ---------------------------------------------------------------------------
void soft_gpu_reinit(const int ntypes, double **cutsq, double **host_prefactor,
                    double **host_cut) {
  int world_me=SLMF.device->world_me();
  int gpu_rank=SLMF.device->gpu_rank();
  int procs_per_gpu=SLMF.device->procs_per_gpu();

  if (world_me==0)
    SLMF.reinit(ntypes, cutsq, host_prefactor, host_cut);

  SLMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      SLMF.reinit(ntypes, cutsq, host_prefactor, host_cut);

    SLMF.device->gpu_barrier();
  }
}

void soft_gpu_clear() {
  SLMF.clear();
}

int ** soft_gpu_compute_n(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success) {
  return SLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                      subhi, tag, nspecial, special, eflag, vflag, eatom,
                      vatom, host_start, ilist, jnum, cpu_time, success);
}

void soft_gpu_compute(const int ago, const int inum_full, const int nall,
                       double **host_x, int *host_type, int *ilist, int *numj,
                       int **firstneigh, const bool eflag, const bool vflag,
                       const bool eatom, const bool vatom, int &host_start,
                       const double cpu_time, bool &success) {
  SLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
               firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double soft_gpu_bytes() {
  return SLMF.host_memory_usage();
}


