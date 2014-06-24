/***************************************************************************
                              coul_debye_ext.cpp
                             -------------------
                              Trung Dac Nguyen

  Functions for LAMMPS access to coul/debye acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : ndtrung@umich.edu
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_coul_debye.h"

using namespace std;
using namespace LAMMPS_AL;

static CoulDebye<PRECISION,ACC_PRECISION> CDEMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int cdebye_gpu_init(const int ntypes, double **host_scale, double **cutsq,
                    double *host_special_coul, const int inum,
                    const int nall, const int max_nbors, const int maxspecial,
                    const double cell_size, int &gpu_mode, FILE *screen,
                    const double qqrd2e, const double kappa) {
  CDEMF.clear();
  gpu_mode=CDEMF.device->gpu_mode();
  double gpu_split=CDEMF.device->particle_split();
  int first_gpu=CDEMF.device->first_device();
  int last_gpu=CDEMF.device->last_device();
  int world_me=CDEMF.device->world_me();
  int gpu_rank=CDEMF.device->gpu_rank();
  int procs_per_gpu=CDEMF.device->procs_per_gpu();

  CDEMF.device->init_message(screen,"coul/debye",first_gpu,last_gpu);

  bool message=false;
  if (CDEMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=CDEMF.init(ntypes, host_scale, cutsq, host_special_coul, inum, nall, 300,
                       maxspecial, cell_size, gpu_split, screen, qqrd2e, kappa);

  CDEMF.device->world_barrier();
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
      init_ok=CDEMF.init(ntypes, host_scale, cutsq, host_special_coul, inum, nall, 300,
                         maxspecial, cell_size, gpu_split, screen, qqrd2e, kappa);

    CDEMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    CDEMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated constants to device
// ---------------------------------------------------------------------------
void cdebye_gpu_reinit(const int ntypes, double **host_scale) {
  int world_me=CDEMF.device->world_me();
  int gpu_rank=CDEMF.device->gpu_rank();
  int procs_per_gpu=CDEMF.device->procs_per_gpu();
  
  if (world_me==0)
    CDEMF.reinit(ntypes, host_scale);
  
  CDEMF.device->world_barrier();
  
  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      CDEMF.reinit(ntypes, host_scale);
    
    CDEMF.device->gpu_barrier();
  }
}

void cdebye_gpu_clear() {
  CDEMF.clear();
}

int** cdebye_gpu_compute_n(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success, double *host_q, double *boxlo,
                           double *prd) {
  return CDEMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, ilist, jnum, cpu_time, success,
                        host_q, boxlo, prd);
}  
			
void cdebye_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success, double *host_q,
                        const int nlocal, double *boxlo, double *prd) {
  CDEMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,firstneigh,eflag,
                vflag,eatom,vatom,host_start,cpu_time,success,host_q,
                nlocal,boxlo,prd);
}

double cdebye_gpu_bytes() {
  return CDEMF.host_memory_usage();
}


