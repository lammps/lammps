/***************************************************************************
                                 born_coul_dsf_ext.cpp
                             -------------------

  Functions for LAMMPS access to born/coul/dsf acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_born_coul_dsf.h"
#include "lammps_gpu.h"

using namespace std;
using namespace LAMMPS_AL;

static BornCoulDSF<PRECISION,ACC_PRECISION> BORNCDMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------

namespace LAMMPS_GPU {
int born_coul_dsf_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                    double **host_born1, double **host_born2, double **host_born3,
                    double **host_a, double **host_c, double **host_d,
                    double **sigma, double **offset, double *special_lj, const int inum,
                    const int nall, const int max_nbors, const int maxspecial,
                    const double cell_size, int &gpu_mode, FILE *screen,
                    double **host_cut_ljsq, double host_cut_coulsq,
                    double *host_special_coul, const double qqrd2e,
                    const double alpha, const double e_shift, const double f_shift) {
  BORNCDMF.clear();
  gpu_mode=BORNCDMF.device->gpu_mode();
  int first_gpu=BORNCDMF.device->first_device();
  int last_gpu=BORNCDMF.device->last_device();
  int world_me=BORNCDMF.device->world_me();
  int gpu_rank=BORNCDMF.device->gpu_rank();
  int procs_per_gpu=BORNCDMF.device->procs_per_gpu();

  BORNCDMF.device->init_message(screen,"born/coul/dsf",first_gpu,last_gpu);

  bool message=false;
  if (BORNCDMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=BORNCDMF.init(ntypes, cutsq, host_rhoinv, host_born1, host_born2,
                          host_born3, host_a, host_c, host_d, sigma,
                          offset, special_lj, inum, nall, max_nbors,
                          maxspecial, cell_size, screen, host_cut_ljsq,
                          host_cut_coulsq, host_special_coul, qqrd2e,
                          alpha, e_shift, f_shift);

  BORNCDMF.device->world_barrier();
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
      init_ok=BORNCDMF.init(ntypes, cutsq, host_rhoinv, host_born1, host_born2,
                            host_born3, host_a, host_c, host_d, sigma,
                            offset, special_lj, inum, nall, max_nbors,
                            maxspecial, cell_size, screen, host_cut_ljsq,
                            host_cut_coulsq, host_special_coul, qqrd2e,
                            alpha, e_shift, f_shift);

    BORNCDMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    BORNCDMF.estimate_gpu_overhead();
  return init_ok;
}

void born_coul_dsf_gpu_clear() {
  BORNCDMF.clear();
}

int **born_coul_dsf_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int **ilist, int **jnum,
                           bool &success, double *host_q, double *boxlo, double *prd,
                           int *periodicity)
{
  return BORNCDMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                          subhi, tag, nspecial, special, eflag, vflag, eatom,
                          vatom, ilist, jnum, success,
                          host_q, boxlo, prd, periodicity);
}

void born_coul_dsf_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, bool &success, double *host_q,
                        const int nlocal, double *boxlo, double *prd) {
  BORNCDMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                   firstneigh,eflag,vflag,eatom,vatom,success,
                   host_q,nlocal,boxlo,prd);
}

double born_coul_dsf_gpu_bytes() {
  return BORNCDMF.host_memory_usage();
}

} // namespace LAMMPS_GPU
