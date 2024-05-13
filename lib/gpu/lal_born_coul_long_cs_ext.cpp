/***************************************************************************
                           born_coul_long_cs_ext.cpp
                             -------------------
                           Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to born/coul/long/cs acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndactrung@gmail.com
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_born_coul_long_cs.h"

using namespace std;
using namespace LAMMPS_AL;

static BornCoulLongCS<PRECISION,ACC_PRECISION> BCLCSMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int bornclcs_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                    double **host_born1, double **host_born2, double **host_born3,
                    double **host_a, double **host_c, double **host_d,
                    double **sigma, double **offset, double *special_lj,
                    const int inum, const int nall, const int max_nbors,
                    const int maxspecial, const double cell_size, int &gpu_mode,
                    FILE *screen, double **host_cut_ljsq, double host_cut_coulsq,
                    double *host_special_coul, const double qqrd2e,
                    const double g_ewald) {
  BCLCSMF.clear();
  gpu_mode=BCLCSMF.device->gpu_mode();
  double gpu_split=BCLCSMF.device->particle_split();
  int first_gpu=BCLCSMF.device->first_device();
  int last_gpu=BCLCSMF.device->last_device();
  int world_me=BCLCSMF.device->world_me();
  int gpu_rank=BCLCSMF.device->gpu_rank();
  int procs_per_gpu=BCLCSMF.device->procs_per_gpu();

  BCLCSMF.device->init_message(screen,"born/coul/long/cs",first_gpu,last_gpu);

  bool message=false;
  if (BCLCSMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=BCLCSMF.init(ntypes, cutsq, host_rhoinv, host_born1, host_born2,
                          host_born3, host_a, host_c, host_d, sigma, offset,
                          special_lj, inum, nall, max_nbors, maxspecial, cell_size,
                          gpu_split, screen, host_cut_ljsq, host_cut_coulsq,
                          host_special_coul, qqrd2e, g_ewald);

  BCLCSMF.device->world_barrier();
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
      init_ok=BCLCSMF.init(ntypes, cutsq, host_rhoinv, host_born1, host_born2,
                            host_born3, host_a, host_c, host_d, sigma, offset,
                            special_lj, inum, nall, max_nbors, maxspecial, cell_size,
                            gpu_split, screen, host_cut_ljsq, host_cut_coulsq,
                            host_special_coul, qqrd2e, g_ewald);

    BCLCSMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    BCLCSMF.estimate_gpu_overhead();
  return init_ok;
}

void bornclcs_gpu_clear() {
  BCLCSMF.clear();
}

int** bornclcs_gpu_compute_n(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum,  const double cpu_time,
                           bool &success, double *host_q, double *boxlo,
                           double *prd) {
  return BCLCSMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                          subhi, tag, nspecial, special, eflag, vflag, eatom,
                          vatom, host_start, ilist, jnum, cpu_time, success,
                          host_q, boxlo, prd);
}

void bornclcs_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success, double *host_q,
                        const int nlocal, double *boxlo, double *prd) {
  BCLCSMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                   firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success,
                   host_q,nlocal,boxlo,prd);
}

double bornclcs_gpu_bytes() {
  return BCLCSMF.host_memory_usage();
}


