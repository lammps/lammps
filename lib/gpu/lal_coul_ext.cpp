/***************************************************************************
                                 coul_ext.cpp
                             -------------------
                               Trung Dac Nguyen

  Functions for LAMMPS access to coul/cut acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndtrung@umich.edu
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_coul.h"

using namespace std;
using namespace LAMMPS_AL;

static Coul<PRECISION,ACC_PRECISION> COULMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int coul_gpu_init(const int ntypes, double **host_scale,
                  double **cutsq, double *special_coul,
                  const int inum, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size,
                  int &gpu_mode, FILE *screen, const double qqrd2e) {
  COULMF.clear();
  gpu_mode=COULMF.device->gpu_mode();
  double gpu_split=COULMF.device->particle_split();
  int first_gpu=COULMF.device->first_device();
  int last_gpu=COULMF.device->last_device();
  int world_me=COULMF.device->world_me();
  int gpu_rank=COULMF.device->gpu_rank();
  int procs_per_gpu=COULMF.device->procs_per_gpu();

  COULMF.device->init_message(screen,"coul/cut",first_gpu,last_gpu);

  bool message=false;
  if (COULMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=COULMF.init(ntypes, host_scale, cutsq, special_coul, inum, nall, 300,
                       maxspecial, cell_size, gpu_split, screen, qqrd2e);

  COULMF.device->world_barrier();
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
      init_ok=COULMF.init(ntypes, host_scale, cutsq, special_coul, inum, nall, 300,
                          maxspecial, cell_size, gpu_split, screen, qqrd2e);

    COULMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    COULMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated constants to device
// ---------------------------------------------------------------------------
void coul_gpu_reinit(const int ntypes, double **host_scale) {
  int world_me=COULMF.device->world_me();
  int gpu_rank=COULMF.device->gpu_rank();
  int procs_per_gpu=COULMF.device->procs_per_gpu();

  if (world_me==0)
    COULMF.reinit(ntypes, host_scale);

  COULMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      COULMF.reinit(ntypes, host_scale);

    COULMF.device->gpu_barrier();
  }
}

void coul_gpu_clear() {
  COULMF.clear();
}

int** coul_gpu_compute_n(const int ago, const int inum_full,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum, const double cpu_time,
                        bool &success, double *host_q, double *boxlo,
                        double *prd) {
  return COULMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success,
                       host_q, boxlo, prd);
}

void coul_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success, double *host_q,
                     const int nlocal, double *boxlo, double *prd) {
  COULMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,firstneigh,eflag,
                vflag,eatom,vatom,host_start,cpu_time,success,host_q,
                nlocal,boxlo,prd);
}

double coul_gpu_bytes() {
  return COULMF.host_memory_usage();
}


