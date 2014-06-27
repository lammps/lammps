/***************************************************************************
                                 gauss_ext.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to gauss acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_gauss.h"

using namespace std;
using namespace LAMMPS_AL;

static Gauss<PRECISION,ACC_PRECISION> GLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int gauss_gpu_init(const int ntypes, double **cutsq, double **host_a, 
                   double **host_b, double **offset, double *special_lj, 
                   const int inum, const int nall, const int max_nbors,  
                   const int maxspecial,
                   const double cell_size, int &gpu_mode, FILE *screen) {
  GLMF.clear();
  gpu_mode=GLMF.device->gpu_mode();
  double gpu_split=GLMF.device->particle_split();
  int first_gpu=GLMF.device->first_device();
  int last_gpu=GLMF.device->last_device();
  int world_me=GLMF.device->world_me();
  int gpu_rank=GLMF.device->gpu_rank();
  int procs_per_gpu=GLMF.device->procs_per_gpu();

  GLMF.device->init_message(screen,"gauss",first_gpu,last_gpu);

  bool message=false;
  if (GLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=GLMF.init(ntypes, cutsq, host_a, host_b, 
                       offset, special_lj, inum, nall, 300,
                       maxspecial, cell_size, gpu_split, screen);

  GLMF.device->world_barrier();
  if (message)
    fprintf(screen,"Done.\n");

  for (int i=0; i<procs_per_gpu; i++) {
    if (message) {
      if (last_gpu-first_gpu==0)
        fprintf(screen,"Initializing GPU %d on core %d...",first_gpu,i);
      else
        fprintf(screen,"Initializing GPUs %d-%d on core %d...",first_gpu,
                last_gpu,i);
      fflush(screen);
    }
    if (gpu_rank==i && world_me!=0)
      init_ok=GLMF.init(ntypes, cutsq, host_a, host_b,
                        offset, special_lj, inum, nall, 300, maxspecial,
                        cell_size, gpu_split, screen);

    GLMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    GLMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void gauss_gpu_reinit(const int ntypes, double **cutsq, double **host_a,
                      double **host_b, double **offset) {
  int world_me=GLMF.device->world_me();
  int gpu_rank=GLMF.device->gpu_rank();
  int procs_per_gpu=GLMF.device->procs_per_gpu();
  
  if (world_me==0)
    GLMF.reinit(ntypes, cutsq, host_a, host_b, offset);
  
  GLMF.device->world_barrier();
  
  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      GLMF.reinit(ntypes, cutsq, host_a, host_b, offset);
    
    GLMF.device->gpu_barrier();
  }
}

void gauss_gpu_clear() {
  GLMF.clear();
}

int ** gauss_gpu_compute_n(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success) {
  return GLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                      subhi, tag, nspecial, special, eflag, vflag, eatom,
                      vatom, host_start, ilist, jnum, cpu_time, success);
}  
			
void gauss_gpu_compute(const int ago, const int inum_full, const int nall,
                       double **host_x, int *host_type, int *ilist, int *numj,
                       int **firstneigh, const bool eflag, const bool vflag,
                       const bool eatom, const bool vatom, int &host_start,
                       const double cpu_time, bool &success) {
  GLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
               firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double gauss_gpu_bytes() {
  return GLMF.host_memory_usage();
}


