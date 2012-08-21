/***************************************************************************
                               dipole_lj_ext.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to dipole/cut acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_dipole_lj.h"

using namespace std;
using namespace LAMMPS_AL;

static DipoleLJ<PRECISION,ACC_PRECISION> DPLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int dpl_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                 double **host_lj2, double **host_lj3, double **host_lj4,
                 double **offset, double *special_lj, const int inum,
                 const int nall, const int max_nbors, const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen,
                 double **host_cut_ljsq, double **host_cut_coulsq,
                 double *host_special_coul, const double qqrd2e) {
  DPLMF.clear();
  gpu_mode=DPLMF.device->gpu_mode();
  double gpu_split=DPLMF.device->particle_split();
  int first_gpu=DPLMF.device->first_device();
  int last_gpu=DPLMF.device->last_device();
  int world_me=DPLMF.device->world_me();
  int gpu_rank=DPLMF.device->gpu_rank();
  int procs_per_gpu=DPLMF.device->procs_per_gpu();

  DPLMF.device->init_message(screen,"dipole/cut",first_gpu,last_gpu);

  bool message=false;
  if (DPLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=DPLMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3,
                       host_lj4, offset, special_lj, inum, nall, 300,
                       maxspecial, cell_size, gpu_split, screen, host_cut_ljsq,
                       host_cut_coulsq, host_special_coul, qqrd2e);

  DPLMF.device->world_barrier();
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
      init_ok=DPLMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                         offset, special_lj, inum, nall, 300, maxspecial,
                         cell_size, gpu_split, screen, host_cut_ljsq,
                         host_cut_coulsq, host_special_coul, qqrd2e);

    DPLMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    DPLMF.estimate_gpu_overhead();
  return init_ok;
}

void dpl_gpu_clear() {
  DPLMF.clear();
}

int** dpl_gpu_compute_n(const int ago, const int inum_full,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, int *tag, int **nspecial, 
                        int **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum, const double cpu_time,
                        bool &success, double *host_q, double **host_mu, 
                        double *boxlo, double *prd) {
  return DPLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success,
                       host_q, host_mu, boxlo, prd);
}  
			
void dpl_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success, double *host_q,
                     double **host_mu, const int nlocal, double *boxlo, double *prd) {
  DPLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,firstneigh,eflag,
                vflag,eatom,vatom,host_start,cpu_time,success,host_q,host_mu,
                nlocal,boxlo,prd);
}

double dpl_gpu_bytes() {
  return DPLMF.host_memory_usage();
}


