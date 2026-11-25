/***************************************************************************
                                 coul_cut_soft_gapsys_ext.cpp
                             -------------------

  Functions for LAMMPS access to coul/cut/soft/gapsys acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_coul_cut_soft_gapsys.h"
#include "lammps_gpu.h"

using namespace std;
using namespace LAMMPS_AL;

static CoulCutSoftGapsys<PRECISION,ACC_PRECISION> CCSGMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------

namespace LAMMPS_GPU {
int coul_cut_soft_gapsys_gpu_init(const int ntypes, double **host_cutsq, double **host_cut_factor,
                                  double *special_lj, double *special_coul, const int inum,
                                  const int nall, const int max_nbors, const int maxspecial,
                                  const double cell_size, int &gpu_mode, FILE *screen,
                                  const double host_qqrd2e, const double host_sigmaq)
{
  CCSGMF.clear();
  gpu_mode=CCSGMF.device->gpu_mode();
  int first_gpu=CCSGMF.device->first_device();
  int last_gpu=CCSGMF.device->last_device();
  int world_me=CCSGMF.device->world_me();
  int gpu_rank=CCSGMF.device->gpu_rank();
  int procs_per_gpu=CCSGMF.device->procs_per_gpu();

  CCSGMF.device->init_message(screen,"coul/cut/soft/gapsys",first_gpu,last_gpu);

  bool message=false;
  if (CCSGMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=CCSGMF.init(ntypes, host_cutsq, host_cut_factor, special_lj, special_coul, inum, nall,
                        max_nbors, maxspecial, cell_size, screen, host_qqrd2e, host_sigmaq);

  CCSGMF.device->world_barrier();
  if (message)
    fprintf(screen,"Done.\n");

  for (int i=0; i<procs_per_gpu; i++) {
    if (message) {
      if (last_gpu-first_gpu==0)
        fprintf(screen,"Initializing Device %d on core %d...",first_gpu,i);
      else
        fprintf(screen,"Initializing Devices %d-%d on core %d...",first_gpu, last_gpu,i);
      fflush(screen);
    }
    if (gpu_rank==i && world_me!=0)
      init_ok=CCSGMF.init(ntypes, host_cutsq, host_cut_factor, special_lj, special_coul, inum, nall,
                          max_nbors, maxspecial, cell_size, screen, host_qqrd2e, host_sigmaq);

    CCSGMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    CCSGMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void coul_cut_soft_gapsys_gpu_reinit(const int ntypes, double **host_cutsq,
                                     double **host_cut_factor) {
  int world_me=CCSGMF.device->world_me();
  int gpu_rank=CCSGMF.device->gpu_rank();
  int procs_per_gpu=CCSGMF.device->procs_per_gpu();

  if (world_me==0)
    CCSGMF.reinit(ntypes, host_cutsq, host_cut_factor);
  CCSGMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      CCSGMF.reinit(ntypes, host_cutsq, host_cut_factor);
    CCSGMF.device->serialize_init();
  }
}

void coul_cut_soft_gapsys_gpu_clear() {
  CCSGMF.clear();
}

int **coul_cut_soft_gapsys_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                         double **host_x, int *host_type, double *sublo,
                                         double *subhi, tagint *tag, int **nspecial,
                                         tagint **special, const bool eflag, const bool vflag,
                                         const bool eatom, const bool vatom, int **ilist,
                                         int **jnum, bool &success, double *host_q, double *boxlo,
                                         double *prd, int *periodicity)
{
  return CCSGMF.compute(ago, inum_full, nall, host_x, host_type, sublo, subhi, tag, nspecial,
                        special, eflag, vflag, eatom, vatom, ilist, jnum, success, host_q, boxlo,
                        prd, periodicity);
}

void coul_cut_soft_gapsys_gpu_compute(const int ago, const int inum_full, const int nall,
                                      double **host_x, int *host_type, int *ilist, int *numj,
                                      int **firstneigh, const bool eflag, const bool vflag,
                                      const bool eatom, const bool vatom, bool &success,
                                      double *host_q, const int nlocal, double *boxlo, double *prd)
{
  CCSGMF.compute(ago, inum_full, nall, host_x, host_type, ilist, numj, firstneigh, eflag, vflag,
                 eatom, vatom, success, host_q, nlocal, boxlo, prd);
}

double coul_cut_soft_gapsys_gpu_bytes() {
  return CCSGMF.host_memory_usage();
}

} // namespace LAMMPS_GPU
