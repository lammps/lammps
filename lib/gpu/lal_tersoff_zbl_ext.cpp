/***************************************************************************
                             tersoff_zbl_ext.cpp
                             -------------------
                              Trung Dac Nguyen

  Functions for LAMMPS access to tersoff/zbl acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndactrung@gmail.com
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_tersoff_zbl.h"

using namespace std;
using namespace LAMMPS_AL;

static TersoffZBL<PRECISION,ACC_PRECISION> TSZMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int tersoff_zbl_gpu_init(const int ntypes, const int inum, const int nall,
                     const int max_nbors, const double cell_size, int &gpu_mode,
                     FILE *screen, int* host_map, const int nelements,
                     int*** host_elem2param, const int nparams,
                     const double* ts_lam1, const double* ts_lam2,
                     const double* ts_lam3, const double* ts_powermint,
                     const double* ts_biga, const double* ts_bigb,
                     const double* ts_bigr, const double* ts_bigd,
                     const double* ts_c1, const double* ts_c2,
                     const double* ts_c3, const double* ts_c4,
                     const double* ts_c, const double* ts_d, const double* ts_h,
                     const double* ts_gamma, const double* ts_beta,
                     const double* ts_powern, const double* ts_Z_i,
                     const double* ts_Z_j, const double* ts_ZBLcut,
                     const double* ts_ZBLexpscale, const double global_e,
                     const double global_a_0, const double global_epsilon_0,
                     const double* ts_cutsq) {
  TSZMF.clear();
  gpu_mode=TSZMF.device->gpu_mode();
  double gpu_split=TSZMF.device->particle_split();
  int first_gpu=TSZMF.device->first_device();
  int last_gpu=TSZMF.device->last_device();
  int world_me=TSZMF.device->world_me();
  int gpu_rank=TSZMF.device->gpu_rank();
  int procs_per_gpu=TSZMF.device->procs_per_gpu();

  // disable host/device split for now
  if (gpu_split != 1.0)
    return -8;

  TSZMF.device->init_message(screen,"tersoff/zbl/gpu",first_gpu,last_gpu);

  bool message=false;
  if (TSZMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=TSZMF.init(ntypes, inum, nall, max_nbors, cell_size, gpu_split, screen,
                      host_map, nelements, host_elem2param, nparams,
                      ts_lam1, ts_lam2, ts_lam3, ts_powermint,
                      ts_biga, ts_bigb, ts_bigr, ts_bigd,
                      ts_c1, ts_c2, ts_c3, ts_c4, ts_c, ts_d, ts_h,
                      ts_gamma, ts_beta, ts_powern, ts_Z_i, ts_Z_j,
                      ts_ZBLcut, ts_ZBLexpscale, global_e, global_a_0,
                      global_epsilon_0, ts_cutsq);

  TSZMF.device->world_barrier();
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
      init_ok=TSZMF.init(ntypes, inum, nall, max_nbors, cell_size, gpu_split, screen,
                        host_map, nelements, host_elem2param, nparams,
                        ts_lam1, ts_lam2, ts_lam3, ts_powermint,
                        ts_biga, ts_bigb, ts_bigr, ts_bigd,
                        ts_c1, ts_c2, ts_c3, ts_c4, ts_c, ts_d, ts_h,
                        ts_gamma, ts_beta, ts_powern, ts_Z_i, ts_Z_j,
                        ts_ZBLcut, ts_ZBLexpscale, global_e, global_a_0,
                        global_epsilon_0, ts_cutsq);

    TSZMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    TSZMF.estimate_gpu_overhead(1);
  return init_ok;
}

void tersoff_zbl_gpu_clear() {
  TSZMF.clear();
}

int ** tersoff_zbl_gpu_compute_n(const int ago, const int inum_full,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum, const double cpu_time,
                        bool &success) {
  return TSZMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success);
}

void tersoff_zbl_gpu_compute(const int ago, const int nlocal, const int nall,
                    const int nlist, double **host_x, int *host_type,
                    int *ilist, int *numj, int **firstneigh, const bool eflag,
                    const bool vflag, const bool eatom, const bool vatom,
                    int &host_start, const double cpu_time, bool &success) {
  TSZMF.compute(ago,nlocal,nall,nlist,host_x,host_type,ilist,numj,
               firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double tersoff_zbl_gpu_bytes() {
  return TSZMF.host_memory_usage();
}


