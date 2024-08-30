/***************************************************************************
                            lj_coul_long_soft_ext.cpp
                            -------------------------
                            Trung Nguyen (U Chicago)

  Functions for LAMMPS access to lj/cut/coul/long/soft acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndactrung@gmail.com
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_lj_coul_long_soft.h"

using namespace std;
using namespace LAMMPS_AL;

static LJCoulLongSoft<PRECISION,ACC_PRECISION> LJCLSMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int ljcls_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                  double **host_lj2, double **host_lj3, double **host_lj4,
                  double **offset, double **epsilon, double *special_lj, const int inum,
                  const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen,
                  double **host_cut_ljsq, double host_cut_coulsq,
                  double *host_special_coul, const double qqrd2e,
                  const double g_ewald) {
  LJCLSMF.clear();
  gpu_mode=LJCLSMF.device->gpu_mode();
  double gpu_split=LJCLSMF.device->particle_split();
  int first_gpu=LJCLSMF.device->first_device();
  int last_gpu=LJCLSMF.device->last_device();
  int world_me=LJCLSMF.device->world_me();
  int gpu_rank=LJCLSMF.device->gpu_rank();
  int procs_per_gpu=LJCLSMF.device->procs_per_gpu();

  LJCLSMF.device->init_message(screen,"lj/cut/coul/long/soft",first_gpu,last_gpu);

  bool message=false;
  if (LJCLSMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=LJCLSMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                        offset, epsilon, special_lj, inum, nall, max_nbors, maxspecial,
                        cell_size, gpu_split, screen, host_cut_ljsq,
                        host_cut_coulsq, host_special_coul, qqrd2e, g_ewald);

  LJCLSMF.device->world_barrier();
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
      init_ok=LJCLSMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                          offset, epsilon, special_lj, inum, nall, max_nbors, maxspecial,
                          cell_size, gpu_split, screen, host_cut_ljsq,
                          host_cut_coulsq, host_special_coul, qqrd2e, g_ewald);

    LJCLSMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    LJCLSMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void ljcls_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1,
                    double **host_lj2, double **host_lj3, double **host_lj4,
                    double **offset, double **epsilon, double **host_cut_ljsq) {
  int world_me=LJCLSMF.device->world_me();
  int gpu_rank=LJCLSMF.device->gpu_rank();
  int procs_per_gpu=LJCLSMF.device->procs_per_gpu();

  if (world_me==0)
    LJCLSMF.reinit(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                  offset, epsilon, host_cut_ljsq);
  LJCLSMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      LJCLSMF.reinit(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                    offset, epsilon, host_cut_ljsq);
    LJCLSMF.device->gpu_barrier();
  }
}

void ljcls_gpu_clear() {
  LJCLSMF.clear();
}

int** ljcls_gpu_compute_n(const int ago, const int inum_full,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum,  const double cpu_time,
                         bool &success, double *host_q, double *boxlo,
                         double *prd) {
  return LJCLSMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, ilist, jnum, cpu_time, success,
                        host_q, boxlo, prd);
}

void ljcls_gpu_compute(const int ago, const int inum_full, const int nall,
                      double **host_x, int *host_type, int *ilist, int *numj,
                      int **firstneigh, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success, double *host_q,
                      const int nlocal, double *boxlo, double *prd) {
  LJCLSMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success,
                host_q,nlocal,boxlo,prd);
}

double ljcls_gpu_bytes() {
  return LJCLSMF.host_memory_usage();
}

