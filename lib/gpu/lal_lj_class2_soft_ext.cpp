/***************************************************************************
                                 lj_class2_soft_ext.cpp
                             -------------------
     Claude Opus 5.5 (Anthropic), under the direction of Axel Kohlmeyer

  Functions for LAMMPS access to lj/class2/soft acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_lj_class2_soft.h"
#include "lammps_gpu.h"

using namespace std;
using namespace LAMMPS_AL;

static LJClass2Soft<PRECISION,ACC_PRECISION> LJC2SOFTMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------

namespace LAMMPS_GPU {
int lj_class2_soft_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                  double **host_lj2, double **host_lj3, double **host_epsilon,
                  double **offset, double *special_lj, const int inum,
                  const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen) {
  LJC2SOFTMF.clear();
  gpu_mode=LJC2SOFTMF.device->gpu_mode();
  int first_gpu=LJC2SOFTMF.device->first_device();
  int last_gpu=LJC2SOFTMF.device->last_device();
  int world_me=LJC2SOFTMF.device->world_me();
  int gpu_rank=LJC2SOFTMF.device->gpu_rank();
  int procs_per_gpu=LJC2SOFTMF.device->procs_per_gpu();

  LJC2SOFTMF.device->init_message(screen,"lj/class2/soft",first_gpu,last_gpu);

  bool message=false;
  if (LJC2SOFTMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=LJC2SOFTMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3,
                        host_epsilon, offset, special_lj, inum, nall, max_nbors,
                        maxspecial, cell_size, screen);

  LJC2SOFTMF.device->world_barrier();
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
      init_ok=LJC2SOFTMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_epsilon,
                          offset, special_lj, inum,  nall, max_nbors, maxspecial,
                          cell_size, screen);

    LJC2SOFTMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    LJC2SOFTMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void lj_class2_soft_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1,
                            double **host_lj2, double **host_lj3, double **host_epsilon,
                            double **offset) {
  int world_me=LJC2SOFTMF.device->world_me();
  int gpu_rank=LJC2SOFTMF.device->gpu_rank();
  int procs_per_gpu=LJC2SOFTMF.device->procs_per_gpu();

  if (world_me==0)
    LJC2SOFTMF.reinit(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_epsilon, offset);
  LJC2SOFTMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      LJC2SOFTMF.reinit(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_epsilon, offset);
    LJC2SOFTMF.device->serialize_init();
  }
}

void lj_class2_soft_gpu_clear() {
  LJC2SOFTMF.clear();
}

int **lj_class2_soft_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum, bool &success, double *prd,
                         int *periodicity)
{
  return LJC2SOFTMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, ilist, jnum, success, prd, periodicity);
}

void lj_class2_soft_gpu_compute(const int ago, const int inum_full, const int nall,
                      double **host_x, int *host_type, int *ilist, int *numj,
                      int **firstneigh, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom, bool &success) {
  LJC2SOFTMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,firstneigh,
                 eflag,vflag,eatom,vatom,success);
}

double lj_class2_soft_gpu_bytes() {
  return LJC2SOFTMF.host_memory_usage();
}

} // namespace LAMMPS_GPU
