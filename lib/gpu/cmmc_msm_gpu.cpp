/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#include <iostream>
#include <cassert>
#include <math.h>

#include "cmmc_msm_gpu_memory.h"

using namespace std;

static CMMM_GPU_Memory<PRECISION,ACC_PRECISION> CMMMMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int cmmm_gpu_init(const int ntypes, double **cutsq, int **cg_type,
                  double **host_lj1, double **host_lj2, double **host_lj3, 
                  double **host_lj4, double **offset, double *special_lj,
                  const int inum, const int nall, const int max_nbors, 
                  const int maxspecial, const double cell_size, int &gpu_mode,
                  FILE *screen, double **host_cut_ljsq, double host_cut_coulsq,
                  double *host_special_coul, const double qqrd2e, 
                  const int smooth) {
  CMMMMF.clear();
  gpu_mode=CMMMMF.device->gpu_mode();
  double gpu_split=CMMMMF.device->particle_split();
  int first_gpu=CMMMMF.device->first_device();
  int last_gpu=CMMMMF.device->last_device();
  int world_me=CMMMMF.device->world_me();
  int gpu_rank=CMMMMF.device->gpu_rank();
  int procs_per_gpu=CMMMMF.device->procs_per_gpu();

  CMMMMF.device->init_message(screen,"cg/cmm/coul/msm",first_gpu,last_gpu);

  bool message=false;
  if (CMMMMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=CMMMMF.init(ntypes, cutsq, cg_type, host_lj1, host_lj2, host_lj3,
                        host_lj4, offset, special_lj, inum,  nall, 300,
                        maxspecial, cell_size, gpu_split, screen, host_cut_ljsq,
                        host_cut_coulsq, host_special_coul, qqrd2e,smooth);

  CMMMMF.device->world_barrier();
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
      init_ok=CMMMMF.init(ntypes, cutsq, cg_type, host_lj1, host_lj2, host_lj3,
                          host_lj4, offset, special_lj, inum,  nall, 300,
                          maxspecial, cell_size, gpu_split, screen,
                          host_cut_ljsq, host_cut_coulsq, host_special_coul,
                          qqrd2e,smooth);

    CMMMMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    CMMMMF.estimate_gpu_overhead();
  return init_ok;
}

void cmmm_gpu_clear() {
  CMMMMF.clear();
}

int** cmmm_gpu_compute_n(const int ago, const int inum_full,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, int *tag, int **nspecial, 
                         int **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum, const double cpu_time,
                         bool &success, double *host_q, double *boxlo,
                         double *prd) {
  return CMMMMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, ilist, jnum, cpu_time, success,
                        host_q, boxlo, prd);
}  
			
void cmmm_gpu_compute(const int ago, const int inum_full, const int nall,
                      double **host_x, int *host_type, int *ilist, int *numj,
                      int **firstneigh, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success, double *host_q,
                      const int nlocal, double *boxlo, double *prd) {
  CMMMMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                 firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success,
                 host_q,nlocal,boxlo,prd);
}

double cmmm_gpu_bytes() {
  return CMMMMF.host_memory_usage();
}


