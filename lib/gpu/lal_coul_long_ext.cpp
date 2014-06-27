/***************************************************************************
                              coul_long_ext.cpp
                             -------------------
                           Axel Kohlmeyer (Temple)

  Functions for LAMMPS access to coul/long acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : July 2011
    email                : a.kohlmeyer@temple.edu
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_coul_long.h"

using namespace std;
using namespace LAMMPS_AL;

static CoulLong<PRECISION,ACC_PRECISION> CLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int cl_gpu_init(const int ntypes, double **host_scale,
                const int inum, const int nall, const int max_nbors,
                const int maxspecial, const double cell_size, int &gpu_mode,
                FILE *screen, double host_cut_coulsq, double *host_special_coul,
                const double qqrd2e, const double g_ewald) {
  CLMF.clear();
  gpu_mode=CLMF.device->gpu_mode();
  double gpu_split=CLMF.device->particle_split();
  int first_gpu=CLMF.device->first_device();
  int last_gpu=CLMF.device->last_device();
  int world_me=CLMF.device->world_me();
  int gpu_rank=CLMF.device->gpu_rank();
  int procs_per_gpu=CLMF.device->procs_per_gpu();

  CLMF.device->init_message(screen,"coul/long",first_gpu,last_gpu);

  bool message=false;
  if (CLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=CLMF.init(ntypes, host_scale, inum, nall, 300, maxspecial,
                      cell_size, gpu_split, screen, host_cut_coulsq,
                      host_special_coul, qqrd2e, g_ewald);

  CLMF.device->world_barrier();
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
      init_ok=CLMF.init(ntypes, host_scale, inum, nall, 300, maxspecial,
                        cell_size, gpu_split, screen, host_cut_coulsq,
                        host_special_coul, qqrd2e, g_ewald);

    CLMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    CLMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void cl_gpu_reinit(const int ntypes, double **host_scale) {
  int world_me=CLMF.device->world_me();
  int gpu_rank=CLMF.device->gpu_rank();
  int procs_per_gpu=CLMF.device->procs_per_gpu();
  
  if (world_me==0)
    CLMF.reinit(ntypes, host_scale);
  
  CLMF.device->world_barrier();
  
  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      CLMF.reinit(ntypes, host_scale);
    
    CLMF.device->gpu_barrier();
  }
}

void cl_gpu_clear() {
  CLMF.clear();
}

int** cl_gpu_compute_n(const int ago, const int inum_full,
		       const int nall, double **host_x, int *host_type,
		       double *sublo, double *subhi, tagint *tag, int **nspecial,
		       tagint **special, const bool eflag, const bool vflag,
		       const bool eatom, const bool vatom, int &host_start,
		       int **ilist, int **jnum,  const double cpu_time,
		       bool &success, double *host_q, double *boxlo,
		       double *prd) {
  return CLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
		      subhi, tag, nspecial, special, eflag, vflag, eatom,
		      vatom, host_start, ilist, jnum, cpu_time, success,
		      host_q, boxlo, prd);
}

void cl_gpu_compute(const int ago, const int inum_full, const int nall,
		    double **host_x, int *host_type, int *ilist, int *numj,
		    int **firstneigh, const bool eflag, const bool vflag,
		    const bool eatom, const bool vatom, int &host_start,
		    const double cpu_time, bool &success, double *host_q,
		    const int nlocal, double *boxlo, double *prd) {
  CLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
	       firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success,
	       host_q,nlocal,boxlo,prd);
}

double cl_gpu_bytes() {
  return CLMF.host_memory_usage();
}


