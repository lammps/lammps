/***************************************************************************
                                 coul_wolf_ext.cpp
                             -------------------

  Functions for LAMMPS access to coul/wolf acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_coul_wolf.h"
#include "lammps_gpu.h"

using namespace std;
using namespace LAMMPS_AL;

static CoulWolf<PRECISION,ACC_PRECISION> CWMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------

namespace LAMMPS_GPU {
int coul_wolf_gpu_init(const int ntypes, const int inum, const int nall,
                  const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen,
                  const double host_cut_coulsq, double *host_special_coul,
                  const double qqrd2e, const double e_shift, const double f_shift,
                  const double alpha) {
  CWMF.clear();
  gpu_mode=CWMF.device->gpu_mode();
  int first_gpu=CWMF.device->first_device();
  int last_gpu=CWMF.device->last_device();
  int world_me=CWMF.device->world_me();
  int gpu_rank=CWMF.device->gpu_rank();
  int procs_per_gpu=CWMF.device->procs_per_gpu();

  CWMF.device->init_message(screen,"coul/wolf",first_gpu,last_gpu);

  bool message=false;
  if (CWMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok = CWMF.init(ntypes, inum, nall, max_nbors, maxspecial, cell_size, screen,
                        host_cut_coulsq, host_special_coul, qqrd2e, e_shift, f_shift, alpha);

  CWMF.device->world_barrier();
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
      init_ok = CWMF.init(ntypes, inum, nall, max_nbors, maxspecial, cell_size, screen,
                          host_cut_coulsq, host_special_coul, qqrd2e, e_shift, f_shift, alpha);

    CWMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    CWMF.estimate_gpu_overhead();
  return init_ok;
}

void coul_wolf_gpu_clear() {
  CWMF.clear();
}

int **coul_wolf_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum, bool &success, double *host_q,
                         double *boxlo, double *prd, int *periodicity)
{
  return CWMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                      subhi, tag, nspecial, special, eflag, vflag, eatom,
                      vatom, ilist, jnum, success,
                      host_q, boxlo, prd, periodicity);
}

void coul_wolf_gpu_compute(const int ago, const int inum_full, const int nall,
                      double **host_x, int *host_type, int *ilist, int *numj,
                      int **firstneigh, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom, bool &success, double *host_q,
                      const int nlocal, double *boxlo, double *prd) {
  CWMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,firstneigh,eflag,
               vflag,eatom,vatom,success,host_q,
               nlocal,boxlo,prd);
}

double coul_wolf_gpu_bytes() {
  return CWMF.host_memory_usage();
}

} // namespace LAMMPS_GPU
