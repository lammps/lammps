/***************************************************************************
                                 lal_dpd_coul_slater_long_ext.cpp
                             -------------------
                            Eddy BARRAUD (IFPEN/Sorbonne)

  Functions for LAMMPS access to dpd/coul/slater/long acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : May 28, 2024
    email                : eddy.barraud@outlook.fr
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_dpd_coul_slater_long.h"

using namespace std;
using namespace LAMMPS_AL;

static DPDCoulSlaterLong<PRECISION,ACC_PRECISION> DPDCMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int dpd_coul_slater_long_gpu_init(const int ntypes, double **host_cutsq, double **host_a0,
                                  double **host_gamma, double **host_sigma, double **host_cut_dpd,
                                  double **host_cut_dpdsq, double **host_cut_slatersq,
                                  double *special_lj, const int inum, const int nall,
                                  const int max_nbors, const int maxspecial,
                                  const double cell_size, int &gpu_mode, FILE *screen,
                                  double *host_special_coul, const double qqrd2e,
                                  const double g_ewald, const double lamda) {
  DPDCMF.clear();
  gpu_mode=DPDCMF.device->gpu_mode();
  double gpu_split=DPDCMF.device->particle_split();
  int first_gpu=DPDCMF.device->first_device();
  int last_gpu=DPDCMF.device->last_device();
  int world_me=DPDCMF.device->world_me();
  int gpu_rank=DPDCMF.device->gpu_rank();
  int procs_per_gpu=DPDCMF.device->procs_per_gpu();

  DPDCMF.device->init_message(screen,"dpd",first_gpu,last_gpu);

  bool message=false;
  if (DPDCMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=DPDCMF.init(ntypes, host_cutsq, host_a0, host_gamma, host_sigma, host_cut_dpd,
                        host_cut_dpdsq, host_cut_slatersq, special_lj, false, inum, nall,
                        max_nbors, maxspecial, cell_size, gpu_split, screen, host_special_coul,
                        qqrd2e, g_ewald, lamda);

  DPDCMF.device->world_barrier();
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
      init_ok=DPDCMF.init(ntypes, host_cutsq, host_a0, host_gamma, host_sigma, host_cut_dpd,
                          host_cut_dpdsq, host_cut_slatersq, special_lj, false, inum, nall,
                          max_nbors, maxspecial, cell_size, gpu_split, screen, host_special_coul,
                          qqrd2e, g_ewald, lamda);

    DPDCMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    DPDCMF.estimate_gpu_overhead();
  return init_ok;
}

void dpd_coul_slater_long_gpu_clear() {
  DPDCMF.clear();
}

int ** dpd_coul_slater_long_gpu_compute_n(const int ago, const int inum_full, const int nall,
                         double **host_x, int *host_type, double *sublo,
                         double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum, const double cpu_time, bool &success,
                         double **host_v, const double dtinvsqrt,
                         const int seed, const int timestep,
                         double *boxlo, double *prd) {
  return DPDCMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success,
                       host_v, dtinvsqrt, seed, timestep, boxlo, prd);
}

void dpd_coul_slater_long_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success, tagint *tag,
                     double **host_v, const double dtinvsqrt,
                     const int seed, const int timestep,
                     const int nlocal, double *boxlo, double *prd) {
  DPDCMF.compute(ago, inum_full, nall, host_x, host_type, ilist, numj,
                firstneigh, eflag, vflag, eatom, vatom, host_start, cpu_time, success,
                tag, host_v, dtinvsqrt, seed, timestep, nlocal, boxlo, prd);
}

void dpd_coul_slater_long_gpu_update_coeff(int ntypes, double **host_a0, double **host_gamma,
                          double **host_sigma, double **host_cut_dpd)
{
   DPDCMF.update_coeff(ntypes,host_a0,host_gamma,host_sigma, host_cut_dpd);
}

void dpd_coul_slater_long_gpu_get_extra_data(double *host_q) {
  DPDCMF.get_extra_data(host_q);
}

double dpd_coul_slater_long_gpu_bytes() {
  return DPDCMF.host_memory_usage();
}


