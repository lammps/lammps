/***************************************************************************
                                 lj_gromacs_coul_gromacs_ext.cpp
                             -------------------

  Functions for LAMMPS access to lj/gromacs/coul/gromacs acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_lj_gromacs_coul_gromacs.h"
#include "lammps_gpu.h"

using namespace std;
using namespace LAMMPS_AL;

static LJGromacsCoulGromacs<PRECISION,ACC_PRECISION> LJGCGMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------

namespace LAMMPS_GPU {
int lj_gromacs_coul_gromacs_gpu_init(const int ntypes, double **host_lj1, double **host_lj2,
                                     double **host_ljsw1, double **host_ljsw2, double **host_lj3,
                                     double **host_lj4, double **host_ljsw3, double **host_ljsw4,
                                     double **host_ljsw5, double *special_lj, double *special_coul,
                                     const int inum, const int nall, const int max_nbors,
                                     const int maxspecial, const double cell_size, int &gpu_mode,
                                     FILE *screen, const double host_qqrd2e,
                                     const double host_cut_bothsq, const double host_cut_ljsq,
                                     const double host_cut_lj_innersq,
                                     const double host_cut_lj_inner, const double host_cut_coulsq,
                                     const double host_cut_coul_innersq,
                                     const double host_cut_coul_inner, const double host_coulsw1,
                                     const double host_coulsw2, const double host_coulsw3,
                                     const double host_coulsw4, const double host_coulsw5)
{
  LJGCGMF.clear();
  gpu_mode=LJGCGMF.device->gpu_mode();
  int first_gpu=LJGCGMF.device->first_device();
  int last_gpu=LJGCGMF.device->last_device();
  int world_me=LJGCGMF.device->world_me();
  int gpu_rank=LJGCGMF.device->gpu_rank();
  int procs_per_gpu=LJGCGMF.device->procs_per_gpu();

  LJGCGMF.device->init_message(screen,"lj/gromacs/coul/gromacs",first_gpu,last_gpu);

  bool message=false;
  if (LJGCGMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=LJGCGMF.init(ntypes, host_lj1, host_lj2, host_ljsw1, host_ljsw2, host_lj3, host_lj4,
                         host_ljsw3, host_ljsw4, host_ljsw5, special_lj, special_coul, inum, nall,
                         max_nbors, maxspecial, cell_size, screen, host_qqrd2e, host_cut_bothsq,
                         host_cut_ljsq, host_cut_lj_innersq, host_cut_lj_inner, host_cut_coulsq,
                         host_cut_coul_innersq, host_cut_coul_inner, host_coulsw1, host_coulsw2,
                         host_coulsw3, host_coulsw4, host_coulsw5);

  LJGCGMF.device->world_barrier();
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
      init_ok=LJGCGMF.init(ntypes, host_lj1, host_lj2, host_ljsw1, host_ljsw2, host_lj3, host_lj4,
                           host_ljsw3, host_ljsw4, host_ljsw5, special_lj, special_coul, inum, nall,
                           max_nbors, maxspecial, cell_size, screen, host_qqrd2e, host_cut_bothsq,
                           host_cut_ljsq, host_cut_lj_innersq, host_cut_lj_inner, host_cut_coulsq,
                           host_cut_coul_innersq, host_cut_coul_inner, host_coulsw1, host_coulsw2,
                           host_coulsw3, host_coulsw4, host_coulsw5);

    LJGCGMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    LJGCGMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void lj_gromacs_coul_gromacs_gpu_reinit(const int ntypes, double **host_lj1, double **host_lj2,
                                        double **host_ljsw1, double **host_ljsw2, double **host_lj3,
                                        double **host_lj4, double **host_ljsw3, double **host_ljsw4,
                                        double **host_ljsw5) {
  int world_me=LJGCGMF.device->world_me();
  int gpu_rank=LJGCGMF.device->gpu_rank();
  int procs_per_gpu=LJGCGMF.device->procs_per_gpu();

  if (world_me==0)
    LJGCGMF.reinit(ntypes, host_lj1, host_lj2, host_ljsw1, host_ljsw2, host_lj3, host_lj4,
                   host_ljsw3, host_ljsw4, host_ljsw5);
  LJGCGMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      LJGCGMF.reinit(ntypes, host_lj1, host_lj2, host_ljsw1, host_ljsw2, host_lj3, host_lj4,
                     host_ljsw3, host_ljsw4, host_ljsw5);
    LJGCGMF.device->serialize_init();
  }
}

void lj_gromacs_coul_gromacs_gpu_clear() {
  LJGCGMF.clear();
}

int **lj_gromacs_coul_gromacs_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                            double **host_x, int *host_type, double *sublo,
                                            double *subhi, tagint *tag, int **nspecial,
                                            tagint **special, const bool eflag, const bool vflag,
                                            const bool eatom, const bool vatom, int **ilist,
                                            int **jnum, bool &success, double *host_q,
                                            double *boxlo, double *prd, int *periodicity)
{
  return LJGCGMF.compute(ago, inum_full, nall, host_x, host_type, sublo, subhi, tag, nspecial,
                         special, eflag, vflag, eatom, vatom, ilist, jnum, success, host_q, boxlo,
                         prd, periodicity);
}

void lj_gromacs_coul_gromacs_gpu_compute(const int ago, const int inum_full, const int nall,
                                         double **host_x, int *host_type, int *ilist, int *numj,
                                         int **firstneigh, const bool eflag, const bool vflag,
                                         const bool eatom, const bool vatom, bool &success,
                                         double *host_q, const int nlocal, double *boxlo,
                                         double *prd)
{
  LJGCGMF.compute(ago, inum_full, nall, host_x, host_type, ilist, numj, firstneigh, eflag, vflag,
                  eatom, vatom, success, host_q, nlocal, boxlo, prd);
}

double lj_gromacs_coul_gromacs_gpu_bytes() {
  return LJGCGMF.host_memory_usage();
}

} // namespace LAMMPS_GPU
