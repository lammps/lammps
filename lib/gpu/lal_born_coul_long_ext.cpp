/***************************************************************************
                           born_coul_long_ext.cpp
                             -------------------
                           Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to born/coul/long acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_born_coul_long.h"

using namespace std;
using namespace LAMMPS_AL;

static BornCoulLong<PRECISION,ACC_PRECISION> BORNCLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int borncl_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                    double **host_born1, double **host_born2, double **host_born3,
                    double **host_a, double **host_c, double **host_d,
                    double **sigma, double **offset, double *special_lj, 
                    const int inum, const int nall, const int max_nbors, 
                    const int maxspecial, const double cell_size, int &gpu_mode, 
                    FILE *screen, double **host_cut_ljsq, double host_cut_coulsq,
                    double *host_special_coul, const double qqrd2e,
                    const double g_ewald) {
  BORNCLMF.clear();
  gpu_mode=BORNCLMF.device->gpu_mode();
  double gpu_split=BORNCLMF.device->particle_split();
  int first_gpu=BORNCLMF.device->first_device();
  int last_gpu=BORNCLMF.device->last_device();
  int world_me=BORNCLMF.device->world_me();
  int gpu_rank=BORNCLMF.device->gpu_rank();
  int procs_per_gpu=BORNCLMF.device->procs_per_gpu();

  BORNCLMF.device->init_message(screen,"born/coul/long",first_gpu,last_gpu);

  bool message=false;
  if (BORNCLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=BORNCLMF.init(ntypes, cutsq, host_rhoinv, host_born1, host_born2, 
                          host_born3, host_a, host_c, host_d, sigma, offset, 
                          special_lj, inum, nall, 300, maxspecial, cell_size, 
                          gpu_split, screen, host_cut_ljsq, host_cut_coulsq, 
                          host_special_coul, qqrd2e, g_ewald);

  BORNCLMF.device->world_barrier();
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
      init_ok=BORNCLMF.init(ntypes, cutsq, host_rhoinv, host_born1, host_born2, 
                            host_born3, host_a, host_c, host_d, sigma, offset, 
                            special_lj, inum, nall, 300, maxspecial, cell_size, 
                            gpu_split, screen, host_cut_ljsq, host_cut_coulsq, 
                            host_special_coul, qqrd2e, g_ewald);

    BORNCLMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    BORNCLMF.estimate_gpu_overhead();
  return init_ok;
}

void borncl_gpu_clear() {
  BORNCLMF.clear();
}

int** borncl_gpu_compute_n(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           double *sublo, double *subhi, int *tag, int **nspecial, 
                           int **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum,  const double cpu_time,
                           bool &success, double *host_q, double *boxlo,
                           double *prd) {
  return BORNCLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                          subhi, tag, nspecial, special, eflag, vflag, eatom,
                          vatom, host_start, ilist, jnum, cpu_time, success,
                          host_q, boxlo, prd);
}  
			
void borncl_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success, double *host_q,
                        const int nlocal, double *boxlo, double *prd) {
  BORNCLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                   firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success,
                   host_q,nlocal,boxlo,prd);
}

double borncl_gpu_bytes() {
  return BORNCLMF.host_memory_usage();
}


