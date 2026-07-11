/***************************************************************************
                             sph_heatconduction_ext.cpp
                             --------------------------
                            Trung Dac Nguyen (U Chicago)

  Functions for LAMMPS access to sph/heatconduction acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : December 2023
    email                : ndactrung@gmail.com
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_sph_heatconduction.h"

using namespace std;
using namespace LAMMPS_AL;

static SPHHeatConduction<PRECISION,ACC_PRECISION> SPHHeatConductionMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int sph_heatconduction_gpu_init(const int ntypes, double **cutsq, double** host_cut,
                    double **host_alpha, double* host_mass, const int dimension,
                    double *special_lj, const int inum, const int nall,
                    const int max_nbors,  const int maxspecial,
                    const double cell_size, int &gpu_mode, FILE *screen) {
  SPHHeatConductionMF.clear();
  gpu_mode=SPHHeatConductionMF.device->gpu_mode();
  double gpu_split=SPHHeatConductionMF.device->particle_split();
  int first_gpu=SPHHeatConductionMF.device->first_device();
  int last_gpu=SPHHeatConductionMF.device->last_device();
  int world_me=SPHHeatConductionMF.device->world_me();
  int gpu_rank=SPHHeatConductionMF.device->gpu_rank();
  int procs_per_gpu=SPHHeatConductionMF.device->procs_per_gpu();

  SPHHeatConductionMF.device->init_message(screen,"sph_heatconduction",first_gpu,last_gpu);

  bool message=false;
  if (SPHHeatConductionMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=SPHHeatConductionMF.init(ntypes, cutsq, host_cut, host_alpha, host_mass,
                         dimension, special_lj, inum, nall, max_nbors,  maxspecial,
                         cell_size, gpu_split, screen);

  SPHHeatConductionMF.device->world_barrier();
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
      init_ok=SPHHeatConductionMF.init(ntypes, cutsq, host_cut, host_alpha, host_mass,
                           dimension, special_lj, inum, nall, max_nbors, maxspecial,
                           cell_size, gpu_split, screen);

    SPHHeatConductionMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    SPHHeatConductionMF.estimate_gpu_overhead();
  return init_ok;
}

void sph_heatconduction_gpu_clear() {
  SPHHeatConductionMF.clear();
}

int ** sph_heatconduction_gpu_compute_n(const int ago, const int inum_full, const int nall,
                            double **host_x, int *host_type, double *sublo,
                            double *subhi, tagint *host_tag, int **nspecial,
                            tagint **special, const bool eflag, const bool vflag,
                            const bool eatom, const bool vatom, int &host_start,
                            int **ilist, int **jnum, const double cpu_time, bool &success,
                            double **host_v) {
  return SPHHeatConductionMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                         subhi, host_tag, nspecial, special, eflag, vflag,
                         eatom, vatom, host_start, ilist, jnum, cpu_time, success,
                         host_v);
}

void sph_heatconduction_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success, tagint *host_tag,
                        double **host_v) {
  SPHHeatConductionMF.compute(ago, inum_full, nall, host_x, host_type, ilist, numj,
                  firstneigh, eflag, vflag, eatom, vatom, host_start, cpu_time, success,
                  host_tag, host_v);
}

void sph_heatconduction_gpu_get_extra_data(double *host_rho, double *host_esph) {
  SPHHeatConductionMF.get_extra_data(host_rho, host_esph);
}

void sph_heatconduction_gpu_update_dE(void **dE_ptr) {
  SPHHeatConductionMF.update_dE(dE_ptr);
}

double sph_heatconduction_gpu_bytes() {
  return SPHHeatConductionMF.host_memory_usage();
}
