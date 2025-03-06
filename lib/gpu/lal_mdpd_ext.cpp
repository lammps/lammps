/***************************************************************************
                                 mdpd_ext.cpp
                             -------------------
                            Trung Dac Nguyen (U Chicago)

  Functions for LAMMPS access to mdpd acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : December 2023
    email                : ndactrung@gmail.com
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_mdpd.h"

using namespace std;
using namespace LAMMPS_AL;

static MDPD<PRECISION,ACC_PRECISION> MDPDMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int mdpd_gpu_init(const int ntypes, double **cutsq,
                  double **host_A_att, double **host_B_rep,
                  double **host_gamma, double **host_sigma,
                  double **host_cut, double **host_cut_r,
                  double *special_lj, const int inum,
                  const int nall, const int max_nbors,  const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen) {
  MDPDMF.clear();
  gpu_mode=MDPDMF.device->gpu_mode();
  double gpu_split=MDPDMF.device->particle_split();
  int first_gpu=MDPDMF.device->first_device();
  int last_gpu=MDPDMF.device->last_device();
  int world_me=MDPDMF.device->world_me();
  int gpu_rank=MDPDMF.device->gpu_rank();
  int procs_per_gpu=MDPDMF.device->procs_per_gpu();

  MDPDMF.device->init_message(screen,"mdpd",first_gpu,last_gpu);

  bool message=false;
  if (MDPDMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=MDPDMF.init(ntypes, cutsq, host_A_att, host_B_rep, host_gamma, host_sigma,
                        host_cut, host_cut_r, special_lj, inum, nall, max_nbors,
                        maxspecial, cell_size, gpu_split, screen);

  MDPDMF.device->world_barrier();
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
      init_ok=MDPDMF.init(ntypes, cutsq, host_A_att, host_B_rep, host_gamma, host_sigma,
                          host_cut, host_cut_r, special_lj, inum, nall, max_nbors,
                          maxspecial, cell_size, gpu_split, screen);

    MDPDMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    MDPDMF.estimate_gpu_overhead();
  return init_ok;
}

void mdpd_gpu_clear() {
  MDPDMF.clear();
}

int ** mdpd_gpu_compute_n(const int ago, const int inum_full, const int nall,
                         double **host_x, int *host_type, double *sublo,
                         double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum, const double cpu_time, bool &success,
                         double **host_v, const double dtinvsqrt,
                         const int seed, const int timestep,
                         double *boxlo, double *prd) {
  return MDPDMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, ilist, jnum, cpu_time, success,
                        host_v, dtinvsqrt, seed, timestep, boxlo, prd);
}

void mdpd_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success, tagint *tag,
                     double **host_v, const double dtinvsqrt,
                     const int seed, const int timestep,
                     const int nlocal, double *boxlo, double *prd) {
  MDPDMF.compute(ago, inum_full, nall, host_x, host_type, ilist, numj,
                 firstneigh, eflag, vflag, eatom, vatom, host_start, cpu_time, success,
                 tag, host_v, dtinvsqrt, seed, timestep, nlocal, boxlo, prd);
}

void mdpd_gpu_get_extra_data(double *host_rho) {
  MDPDMF.get_extra_data(host_rho);
}

double mdpd_gpu_bytes() {
  return MDPDMF.host_memory_usage();
}


