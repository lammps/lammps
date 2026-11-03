/***************************************************************************
                                 lj_coul_soft_gapsys_ext.cpp
                             -------------------

  Functions for LAMMPS access to lj/cut/coul/cut/soft/gapsys acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_lj_coul_soft_gapsys.h"
#include "lammps_gpu.h"

using namespace std;
using namespace LAMMPS_AL;

static LJCoulSoftGapsys<PRECISION,ACC_PRECISION> LJCCSGMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------

namespace LAMMPS_GPU {
int lj_cut_coul_cut_soft_gapsys_gpu_init(const int ntypes, double **host_lj1, double **host_lj2,
                                         double **host_cut_ljsq, double **host_cut_lj_inner,
                                         double **host_coeff_b1, double **host_coeff_b2,
                                         double **host_cut_coulsq, double **host_cut_factor,
                                         double **host_cutsq, double **host_lj3, double **host_lj4,
                                         double **host_offset, double **host_coeff_a3,
                                         double *special_lj, double *special_coul, const int inum,
                                         const int nall, const int max_nbors, const int maxspecial,
                                         const double cell_size, int &gpu_mode, FILE *screen,
                                         const double host_qqrd2e, const double host_sigmaq)
{
  LJCCSGMF.clear();
  gpu_mode=LJCCSGMF.device->gpu_mode();
  int first_gpu=LJCCSGMF.device->first_device();
  int last_gpu=LJCCSGMF.device->last_device();
  int world_me=LJCCSGMF.device->world_me();
  int gpu_rank=LJCCSGMF.device->gpu_rank();
  int procs_per_gpu=LJCCSGMF.device->procs_per_gpu();

  LJCCSGMF.device->init_message(screen,"lj/cut/coul/cut/soft/gapsys",first_gpu,last_gpu);

  bool message=false;
  if (LJCCSGMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=LJCCSGMF.init(ntypes, host_lj1, host_lj2, host_cut_ljsq, host_cut_lj_inner,
                          host_coeff_b1, host_coeff_b2, host_cut_coulsq, host_cut_factor,
                          host_cutsq, host_lj3, host_lj4, host_offset, host_coeff_a3, special_lj,
                          special_coul, inum, nall, max_nbors, maxspecial, cell_size, screen,
                          host_qqrd2e, host_sigmaq);

  LJCCSGMF.device->world_barrier();
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
      init_ok=LJCCSGMF.init(ntypes, host_lj1, host_lj2, host_cut_ljsq, host_cut_lj_inner,
                            host_coeff_b1, host_coeff_b2, host_cut_coulsq, host_cut_factor,
                            host_cutsq, host_lj3, host_lj4, host_offset, host_coeff_a3, special_lj,
                            special_coul, inum, nall, max_nbors, maxspecial, cell_size, screen,
                            host_qqrd2e, host_sigmaq);

    LJCCSGMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    LJCCSGMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void lj_cut_coul_cut_soft_gapsys_gpu_reinit(const int ntypes, double **host_lj1, double **host_lj2,
                                            double **host_cut_ljsq, double **host_cut_lj_inner,
                                            double **host_coeff_b1, double **host_coeff_b2,
                                            double **host_cut_coulsq, double **host_cut_factor,
                                            double **host_cutsq, double **host_lj3,
                                            double **host_lj4, double **host_offset,
                                            double **host_coeff_a3) {
  int world_me=LJCCSGMF.device->world_me();
  int gpu_rank=LJCCSGMF.device->gpu_rank();
  int procs_per_gpu=LJCCSGMF.device->procs_per_gpu();

  if (world_me==0)
    LJCCSGMF.reinit(ntypes, host_lj1, host_lj2, host_cut_ljsq, host_cut_lj_inner, host_coeff_b1,
                    host_coeff_b2, host_cut_coulsq, host_cut_factor, host_cutsq, host_lj3, host_lj4,
                    host_offset, host_coeff_a3);
  LJCCSGMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      LJCCSGMF.reinit(ntypes, host_lj1, host_lj2, host_cut_ljsq, host_cut_lj_inner, host_coeff_b1,
                      host_coeff_b2, host_cut_coulsq, host_cut_factor, host_cutsq, host_lj3,
                      host_lj4, host_offset, host_coeff_a3);
    LJCCSGMF.device->serialize_init();
  }
}

void lj_cut_coul_cut_soft_gapsys_gpu_clear() {
  LJCCSGMF.clear();
}

int **lj_cut_coul_cut_soft_gapsys_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                                double **host_x, int *host_type, double *sublo,
                                                double *subhi, tagint *tag, int **nspecial,
                                                tagint **special, const bool eflag,
                                                const bool vflag, const bool eatom,
                                                const bool vatom, int **ilist, int **jnum,
                                                bool &success, double *host_q, double *boxlo,
                                                double *prd, int *periodicity)
{
  return LJCCSGMF.compute(ago, inum_full, nall, host_x, host_type, sublo, subhi, tag, nspecial,
                          special, eflag, vflag, eatom, vatom, ilist, jnum, success, host_q, boxlo,
                          prd, periodicity);
}

void lj_cut_coul_cut_soft_gapsys_gpu_compute(const int ago, const int inum_full, const int nall,
                                             double **host_x, int *host_type, int *ilist, int *numj,
                                             int **firstneigh, const bool eflag, const bool vflag,
                                             const bool eatom, const bool vatom, bool &success,
                                             double *host_q, const int nlocal, double *boxlo,
                                             double *prd)
{
  LJCCSGMF.compute(ago, inum_full, nall, host_x, host_type, ilist, numj, firstneigh, eflag, vflag,
                   eatom, vatom, success, host_q, nlocal, boxlo, prd);
}

double lj_cut_coul_cut_soft_gapsys_gpu_bytes() {
  return LJCCSGMF.host_memory_usage();
}

} // namespace LAMMPS_GPU
