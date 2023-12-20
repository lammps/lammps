/***************************************************************************
                                 edpd_ext.cpp
                             -------------------
                            Trung Dac Nguyen (U Chicago)

  Functions for LAMMPS access to edpd acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : September 2023
    email                : ndactrung@gmail.com
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_edpd.h"

using namespace std;
using namespace LAMMPS_AL;

static EDPD<PRECISION,ACC_PRECISION> EDPDMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int edpd_gpu_init(const int ntypes, double **cutsq, double **host_a0,
                  double **host_gamma, double **host_cut, double **host_power,
                  double **host_kappa, double **host_powerT, double **host_cutT,
                  double ***host_sc, double ***host_kc, double *host_mass,
                  double *special_lj, const int power_flag, const int kappa_flag,
                  const int inum, const int nall,
                  const int max_nbors,  const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen) {
  EDPDMF.clear();
  gpu_mode=EDPDMF.device->gpu_mode();
  double gpu_split=EDPDMF.device->particle_split();
  int first_gpu=EDPDMF.device->first_device();
  int last_gpu=EDPDMF.device->last_device();
  int world_me=EDPDMF.device->world_me();
  int gpu_rank=EDPDMF.device->gpu_rank();
  int procs_per_gpu=EDPDMF.device->procs_per_gpu();

  EDPDMF.device->init_message(screen,"edpd",first_gpu,last_gpu);

  bool message=false;
  if (EDPDMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=EDPDMF.init(ntypes, cutsq, host_a0, host_gamma, host_cut,
                        host_power, host_kappa, host_powerT,
                        host_cutT, host_sc, host_kc, host_mass,
                        special_lj, power_flag, kappa_flag,
                        inum, nall, max_nbors,  maxspecial,
                        cell_size, gpu_split, screen);

  EDPDMF.device->world_barrier();
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
      init_ok=EDPDMF.init(ntypes, cutsq, host_a0, host_gamma, host_cut,
                          host_power, host_kappa, host_powerT, host_cutT,
                          host_sc, host_kc, host_mass,
                          special_lj, power_flag, kappa_flag,
                          inum, nall, max_nbors, maxspecial,
                          cell_size, gpu_split, screen);

    EDPDMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    EDPDMF.estimate_gpu_overhead();
  return init_ok;
}

void edpd_gpu_clear() {
  EDPDMF.clear();
}

int ** edpd_gpu_compute_n(const int ago, const int inum_full, const int nall,
                         double **host_x, int *host_type, double *sublo,
                         double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum, const double cpu_time, bool &success,
                         double **host_v, const double dtinvsqrt,
                         const int seed, const int timestep,
                         double *boxlo, double *prd) {
  return EDPDMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, ilist, jnum, cpu_time, success,
                        host_v, dtinvsqrt, seed, timestep, boxlo, prd);
}

void edpd_gpu_compute(const int ago, const int inum_full, const int nall,
                      double **host_x, int *host_type, int *ilist, int *numj,
                      int **firstneigh, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success, tagint *tag,
                      double **host_v, const double dtinvsqrt,
                      const int seed, const int timestep,
                      const int nlocal, double *boxlo, double *prd) {
  EDPDMF.compute(ago, inum_full, nall, host_x, host_type, ilist, numj,
                firstneigh, eflag, vflag, eatom, vatom, host_start, cpu_time, success,
                tag, host_v, dtinvsqrt, seed, timestep, nlocal, boxlo, prd);
}

void edpd_gpu_get_extra_data(double *host_T, double *host_cv) {
  EDPDMF.get_extra_data(host_T, host_cv);
}

void edpd_gpu_update_flux(void **flux_ptr) {
  EDPDMF.update_flux(flux_ptr);
}

double edpd_gpu_bytes() {
  return EDPDMF.host_memory_usage();
}
