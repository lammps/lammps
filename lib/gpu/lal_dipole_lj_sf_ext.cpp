/***************************************************************************
                               dipole_lj_sf_ext.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to dipole/sf acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_dipole_lj_sf.h"

using namespace std;
using namespace LAMMPS_AL;

static DipoleLJSF<PRECISION,ACC_PRECISION> DPLSFMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int dplsf_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                   double **host_lj2, double **host_lj3, double **host_lj4,
                   double *special_lj, const int inum,
                   const int nall, const int max_nbors, const int maxspecial,
                   const double cell_size, int &gpu_mode, FILE *screen,
                   double **host_cut_ljsq, double **host_cut_coulsq,
                   double *host_special_coul, const double qqrd2e) {
  DPLSFMF.clear();
  gpu_mode=DPLSFMF.device->gpu_mode();
  double gpu_split=DPLSFMF.device->particle_split();
  int first_gpu=DPLSFMF.device->first_device();
  int last_gpu=DPLSFMF.device->last_device();
  int world_me=DPLSFMF.device->world_me();
  int gpu_rank=DPLSFMF.device->gpu_rank();
  int procs_per_gpu=DPLSFMF.device->procs_per_gpu();

  DPLSFMF.device->init_message(screen,"dipole/sf",first_gpu,last_gpu);

  bool message=false;
  if (DPLSFMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=DPLSFMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3,
                         host_lj4, special_lj, inum, nall, 300,
                         maxspecial, cell_size, gpu_split, screen, host_cut_ljsq,
                         host_cut_coulsq, host_special_coul, qqrd2e);

  DPLSFMF.device->world_barrier();
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
      init_ok=DPLSFMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                           special_lj, inum, nall, 300, maxspecial,
                           cell_size, gpu_split, screen, host_cut_ljsq,
                           host_cut_coulsq, host_special_coul, qqrd2e);

    DPLSFMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    DPLSFMF.estimate_gpu_overhead();
  return init_ok;
}

void dplsf_gpu_clear() {
  DPLSFMF.clear();
}

int** dplsf_gpu_compute_n(const int ago, const int inum_full,
                          const int nall, double **host_x, int *host_type,
                          double *sublo, double *subhi, tagint *tag, int **nspecial, 
                          tagint **special, const bool eflag, const bool vflag,
                          const bool eatom, const bool vatom, int &host_start,
                          int **ilist, int **jnum, const double cpu_time,
                          bool &success, double *host_q, double **host_mu, 
                          double *boxlo, double *prd) {
  return DPLSFMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                         subhi, tag, nspecial, special, eflag, vflag, eatom,
                         vatom, host_start, ilist, jnum, cpu_time, success,
                         host_q, host_mu, boxlo, prd);
}  
			
void dplsf_gpu_compute(const int ago, const int inum_full, const int nall,
                       double **host_x, int *host_type, int *ilist, int *numj,
                       int **firstneigh, const bool eflag, const bool vflag,
                       const bool eatom, const bool vatom, int &host_start,
                       const double cpu_time, bool &success, double *host_q,
                       double **host_mu, const int nlocal, double *boxlo, double *prd) {
  DPLSFMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,firstneigh,eflag,
                  vflag,eatom,vatom,host_start,cpu_time,success,host_q,host_mu,
                  nlocal,boxlo,prd);
}

double dplsf_gpu_bytes() {
  return DPLSFMF.host_memory_usage();
}


