/***************************************************************************
                               lj_cubic_ext.cpp
                             -------------------
                               Trung Dac Nguyen

  Functions for LAMMPS access to lj/cubic acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndactrung@gmail.com
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_lj_cubic.h"

using namespace std;
using namespace LAMMPS_AL;

static LJCubic<PRECISION,ACC_PRECISION> LJCubicLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int ljcb_gpu_init(const int ntypes, double **cutsq, double **cut_inner_sq,
                  double **cut_inner, double **sigma, double **epsilon,
                  double **host_lj1, double **host_lj2, double **host_lj3,
                  double **host_lj4, double *special_lj,
                  const int inum, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size,
                  int &gpu_mode, FILE *screen) {
  LJCubicLMF.clear();
  gpu_mode=LJCubicLMF.device->gpu_mode();
  double gpu_split=LJCubicLMF.device->particle_split();
  int first_gpu=LJCubicLMF.device->first_device();
  int last_gpu=LJCubicLMF.device->last_device();
  int world_me=LJCubicLMF.device->world_me();
  int gpu_rank=LJCubicLMF.device->gpu_rank();
  int procs_per_gpu=LJCubicLMF.device->procs_per_gpu();

  LJCubicLMF.device->init_message(screen,"lj/cubic",first_gpu,last_gpu);

  bool message=false;
  if (LJCubicLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=LJCubicLMF.init(ntypes, cutsq, cut_inner_sq, cut_inner, sigma,
                            epsilon, host_lj1, host_lj2, host_lj3, host_lj4,
                            special_lj, inum, nall, 300, maxspecial,
                            cell_size, gpu_split, screen);

  LJCubicLMF.device->world_barrier();
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
      init_ok=LJCubicLMF.init(ntypes, cutsq, cut_inner_sq, cut_inner, sigma,
                              epsilon, host_lj1, host_lj2, host_lj3, host_lj4,
                              special_lj, inum, nall, 300, maxspecial,
                              cell_size, gpu_split, screen);

    LJCubicLMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    LJCubicLMF.estimate_gpu_overhead();
  return init_ok;
}

void ljcb_gpu_clear() {
  LJCubicLMF.clear();
}

int ** ljcb_gpu_compute_n(const int ago, const int inum_full,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum, const double cpu_time,
                         bool &success) {
  return LJCubicLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success);
}

void ljcb_gpu_compute(const int ago, const int inum_full, const int nall,
                      double **host_x, int *host_type, int *ilist, int *numj,
                      int **firstneigh, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success) {
  LJCubicLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double ljcb_gpu_bytes() {
  return LJCubicLMF.host_memory_usage();
}


