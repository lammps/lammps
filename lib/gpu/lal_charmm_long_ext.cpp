/***************************************************************************
                             charmm_long_ext.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Functions for LAMMPS access to charmm/coul/long acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_charmm_long.h"

using namespace std;
using namespace LAMMPS_AL;

static CHARMMLong<PRECISION,ACC_PRECISION> CRMLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int crml_gpu_init(const int ntypes, double cut_bothsq, double **host_lj1,
                  double **host_lj2, double **host_lj3, double **host_lj4,
                  double **offset, double *special_lj, const int inum,
                  const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen,
                  double host_cut_ljsq, double host_cut_coulsq,
                  double *host_special_coul, const double qqrd2e,
                  const double g_ewald, const double cut_lj_innersq,
                  const double denom_lj, double **epsilon,
                  double **sigma, const bool mix_arithmetic) {
  CRMLMF.clear();
  gpu_mode=CRMLMF.device->gpu_mode();
  double gpu_split=CRMLMF.device->particle_split();
  int first_gpu=CRMLMF.device->first_device();
  int last_gpu=CRMLMF.device->last_device();
  int world_me=CRMLMF.device->world_me();
  int gpu_rank=CRMLMF.device->gpu_rank();
  int procs_per_gpu=CRMLMF.device->procs_per_gpu();

  CRMLMF.device->init_message(screen,"lj/charmm/coul/long",first_gpu,last_gpu);

  bool message=false;
  if (CRMLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    CRMLMF.init(ntypes, cut_bothsq, host_lj1, host_lj2, host_lj3, host_lj4,
                offset, special_lj, inum, nall, 300, maxspecial, cell_size,
                gpu_split, screen, host_cut_ljsq, host_cut_coulsq,
                host_special_coul, qqrd2e, g_ewald, cut_lj_innersq, denom_lj,
                epsilon,sigma,mix_arithmetic);

  CRMLMF.device->world_barrier();
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
      init_ok=CRMLMF.init(ntypes, cut_bothsq, host_lj1, host_lj2, host_lj3,
                          host_lj4, offset, special_lj, inum, nall, 300,
                          maxspecial, cell_size, gpu_split, screen,
                          host_cut_ljsq, host_cut_coulsq, host_special_coul,
                          qqrd2e, g_ewald,  cut_lj_innersq, denom_lj, epsilon,
                          sigma, mix_arithmetic);

    CRMLMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    CRMLMF.estimate_gpu_overhead();
  return init_ok;
}

void crml_gpu_clear() {
  CRMLMF.clear();
}

int** crml_gpu_compute_n(const int ago, const int inum_full,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, tagint *tag, int **nspecial, 
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum, const double cpu_time,
                         bool &success, double *host_q, double *boxlo,
                         double *prd) {
  return CRMLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, ilist, jnum, cpu_time, success,
                        host_q, boxlo, prd);
}  
			
void crml_gpu_compute(const int ago, const int inum_full,
	 	                  const int nall, double **host_x, int *host_type,
                      int *ilist, int *numj, int **firstneigh,
		                  const bool eflag, const bool vflag, const bool eatom,
                      const bool vatom, int &host_start, const double cpu_time,
                      bool &success, double *host_q, const int nlocal, 
                      double *boxlo, double *prd) {
  CRMLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,firstneigh,
                 eflag,vflag,eatom,vatom,host_start,cpu_time,success,host_q,
                 nlocal,boxlo,prd);
}

double crml_gpu_bytes() {
  return CRMLMF.host_memory_usage();
}


