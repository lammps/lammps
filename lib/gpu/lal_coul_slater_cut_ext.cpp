/***************************************************************************
                                 coul_slater_cut_ext.cpp
                             -------------------
     Claude Opus 5.5 (Anthropic), under the direction of Axel Kohlmeyer

  Functions for LAMMPS access to coul/slater/cut acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_coul_slater_cut.h"
#include "lammps_gpu.h"

using namespace std;
using namespace LAMMPS_AL;

static CoulSlaterCut<PRECISION,ACC_PRECISION> CSCMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------

namespace LAMMPS_GPU {
int coul_slater_cut_gpu_init(const int ntypes, double **host_scale,
                  double **cutsq, double *special_coul,
                  const int inum, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size,
                  int &gpu_mode, FILE *screen, const double qqrd2e,
                  const double lamda) {
  CSCMF.clear();
  gpu_mode=CSCMF.device->gpu_mode();
  int first_gpu=CSCMF.device->first_device();
  int last_gpu=CSCMF.device->last_device();
  int world_me=CSCMF.device->world_me();
  int gpu_rank=CSCMF.device->gpu_rank();
  int procs_per_gpu=CSCMF.device->procs_per_gpu();

  CSCMF.device->init_message(screen,"coul/cut",first_gpu,last_gpu);

  bool message=false;
  if (CSCMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=CSCMF.init(ntypes, host_scale, cutsq, special_coul, inum, nall, max_nbors,
                       maxspecial, cell_size, screen, qqrd2e, lamda);

  CSCMF.device->world_barrier();
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
      init_ok=CSCMF.init(ntypes, host_scale, cutsq, special_coul, inum, nall, max_nbors,
                          maxspecial, cell_size, screen, qqrd2e, lamda);

    CSCMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    CSCMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated constants to device
// ---------------------------------------------------------------------------
void coul_slater_cut_gpu_reinit(const int ntypes, double **host_scale) {
  int world_me=CSCMF.device->world_me();
  int gpu_rank=CSCMF.device->gpu_rank();
  int procs_per_gpu=CSCMF.device->procs_per_gpu();

  if (world_me==0)
    CSCMF.reinit(ntypes, host_scale);

  CSCMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      CSCMF.reinit(ntypes, host_scale);

    CSCMF.device->serialize_init();
  }
}

void coul_slater_cut_gpu_clear() {
  CSCMF.clear();
}

int **coul_slater_cut_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum, bool &success, double *host_q,
                         double *boxlo, double *prd, int *periodicity)
{
  return CSCMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, ilist, jnum, success,
                       host_q, boxlo, prd, periodicity);
}

void coul_slater_cut_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, bool &success, double *host_q,
                     const int nlocal, double *boxlo, double *prd) {
  CSCMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,firstneigh,eflag,
                vflag,eatom,vatom,success,host_q,
                nlocal,boxlo,prd);
}

double coul_slater_cut_gpu_bytes() {
  return CSCMF.host_memory_usage();
}

} // namespace LAMMPS_GPU
