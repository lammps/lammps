/***************************************************************************
                                 coul_long_soft_ext.cpp
                             -------------------
     Claude Opus 5.5 (Anthropic), under the direction of Axel Kohlmeyer

  Functions for LAMMPS access to coul/long/soft acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_coul_long_soft.h"
#include "lammps_gpu.h"

using namespace std;
using namespace LAMMPS_AL;

static CoulLongSoft<PRECISION,ACC_PRECISION> CLSMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------

namespace LAMMPS_GPU {
int coul_long_soft_gpu_init(const int ntypes, double **host_scale, double **host_lam1,
                double **host_lam2,
                const int inum, const int nall, const int max_nbors,
                const int maxspecial, const double cell_size, int &gpu_mode,
                FILE *screen, double host_cut_coulsq, double *host_special_coul,
                const double qqrd2e, const double g_ewald) {
  CLSMF.clear();
  gpu_mode=CLSMF.device->gpu_mode();
  int first_gpu=CLSMF.device->first_device();
  int last_gpu=CLSMF.device->last_device();
  int world_me=CLSMF.device->world_me();
  int gpu_rank=CLSMF.device->gpu_rank();
  int procs_per_gpu=CLSMF.device->procs_per_gpu();

  CLSMF.device->init_message(screen,"coul/long/soft",first_gpu,last_gpu);

  bool message=false;
  if (CLSMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=CLSMF.init(ntypes, host_scale, host_lam1, host_lam2, inum, nall, max_nbors, maxspecial,
                      cell_size, screen, host_cut_coulsq,
                      host_special_coul, qqrd2e, g_ewald);

  CLSMF.device->world_barrier();
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
      init_ok=CLSMF.init(ntypes, host_scale, host_lam1, host_lam2, inum, nall, max_nbors, maxspecial,
                        cell_size, screen, host_cut_coulsq,
                        host_special_coul, qqrd2e, g_ewald);

    CLSMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    CLSMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void coul_long_soft_gpu_reinit(const int ntypes, double **host_scale, double **host_lam1,
                             double **host_lam2) {
  int world_me=CLSMF.device->world_me();
  int gpu_rank=CLSMF.device->gpu_rank();
  int procs_per_gpu=CLSMF.device->procs_per_gpu();

  if (world_me==0)
    CLSMF.reinit(ntypes, host_scale, host_lam1, host_lam2);

  CLSMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      CLSMF.reinit(ntypes, host_scale, host_lam1, host_lam2);

    CLSMF.device->serialize_init();
  }
}

void coul_long_soft_gpu_clear() {
  CLSMF.clear();
}

int **coul_long_soft_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                       int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                       tagint **special, const bool eflag, const bool vflag, const bool eatom,
                       const bool vatom, int **ilist, int **jnum, bool &success, double *host_q,
                       double *boxlo, double *prd, int *periodicity)
{
  return CLSMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                      subhi, tag, nspecial, special, eflag, vflag, eatom,
                      vatom, ilist, jnum, success,
                      host_q, boxlo, prd, periodicity);
}

void coul_long_soft_gpu_compute(const int ago, const int inum_full, const int nall,
                    double **host_x, int *host_type, int *ilist, int *numj,
                    int **firstneigh, const bool eflag, const bool vflag,
                    const bool eatom, const bool vatom, bool &success, double *host_q,
                    const int nlocal, double *boxlo, double *prd) {
  CLSMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
               firstneigh,eflag,vflag,eatom,vatom,success,
               host_q,nlocal,boxlo,prd);
}

double coul_long_soft_gpu_bytes() {
  return CLSMF.host_memory_usage();
}

} // namespace LAMMPS_GPU
