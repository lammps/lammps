/***************************************************************************
                           born_coul_wolf_cs_ext.cpp
                             -------------------
                           Trung Dac Nguyen (Northwestern)

  Functions for LAMMPS access to born/coul/wolf/cs acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndactrung@gmail.com
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_born_coul_wolf_cs.h"

using namespace std;
using namespace LAMMPS_AL;

static BornCoulWolfCS<PRECISION,ACC_PRECISION> BornCWCST;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int borncwcs_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                    double **host_born1, double **host_born2, double **host_born3,
                    double **host_a, double **host_c, double **host_d,
                    double **sigma, double **offset, double *special_lj, const int inum,
                    const int nall, const int max_nbors, const int maxspecial,
                    const double cell_size, int &gpu_mode, FILE *screen,
                    double **host_cut_ljsq, double host_cut_coulsq,
                    double *host_special_coul, const double qqrd2e,
                    const double alf, const double e_shift, const double f_shift) {
  BornCWCST.clear();
  gpu_mode=BornCWCST.device->gpu_mode();
  double gpu_split=BornCWCST.device->particle_split();
  int first_gpu=BornCWCST.device->first_device();
  int last_gpu=BornCWCST.device->last_device();
  int world_me=BornCWCST.device->world_me();
  int gpu_rank=BornCWCST.device->gpu_rank();
  int procs_per_gpu=BornCWCST.device->procs_per_gpu();

  BornCWCST.device->init_message(screen,"born/coul/wolf/cs",first_gpu,last_gpu);

  bool message=false;
  if (BornCWCST.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=BornCWCST.init(ntypes, cutsq, host_rhoinv, host_born1, host_born2,
                          host_born3, host_a, host_c, host_d, sigma,
                          offset, special_lj, inum, nall, 300,
                          maxspecial, cell_size, gpu_split, screen, host_cut_ljsq,
                          host_cut_coulsq, host_special_coul, qqrd2e,
                          alf, e_shift, f_shift);

  BornCWCST.device->world_barrier();
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
      init_ok=BornCWCST.init(ntypes, cutsq, host_rhoinv, host_born1, host_born2,
                            host_born3, host_a, host_c, host_d, sigma,
                            offset, special_lj, inum, nall, 300,
                            maxspecial, cell_size, gpu_split, screen, host_cut_ljsq,
                            host_cut_coulsq, host_special_coul, qqrd2e,
                            alf, e_shift, f_shift);

    BornCWCST.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    BornCWCST.estimate_gpu_overhead();
  return init_ok;
}

void borncwcs_gpu_clear() {
  BornCWCST.clear();
}

int** borncwcs_gpu_compute_n(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum,  const double cpu_time,
                           bool &success, double *host_q, double *boxlo,
                           double *prd) {
  return BornCWCST.compute(ago, inum_full, nall, host_x, host_type, sublo,
                          subhi, tag, nspecial, special, eflag, vflag, eatom,
                          vatom, host_start, ilist, jnum, cpu_time, success,
                          host_q, boxlo, prd);
}

void borncwcs_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success, double *host_q,
                        const int nlocal, double *boxlo, double *prd) {
  BornCWCST.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                   firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success,
                   host_q,nlocal,boxlo,prd);
}

double borncwcs_gpu_bytes() {
  return BornCWCST.host_memory_usage();
}


