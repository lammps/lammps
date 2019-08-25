/***************************************************************************
                                 dpd_ext.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to dpd acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Jan 15, 2014
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_dpd.h"

using namespace std;
using namespace LAMMPS_AL;

static DPD<PRECISION,ACC_PRECISION> DPDMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int dpd_gpu_init(const int ntypes, double **cutsq, double **host_a0,
                 double **host_gamma, double **host_sigma, double **host_cut,
                 double *special_lj, bool tstat_only, const int inum,
                 const int nall, const int max_nbors,  const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen) {
  DPDMF.clear();
  gpu_mode=DPDMF.device->gpu_mode();
  double gpu_split=DPDMF.device->particle_split();
  int first_gpu=DPDMF.device->first_device();
  int last_gpu=DPDMF.device->last_device();
  int world_me=DPDMF.device->world_me();
  int gpu_rank=DPDMF.device->gpu_rank();
  int procs_per_gpu=DPDMF.device->procs_per_gpu();

  DPDMF.device->init_message(screen,"dpd",first_gpu,last_gpu);

  bool message=false;
  if (DPDMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=DPDMF.init(ntypes, cutsq, host_a0, host_gamma, host_sigma,
                       host_cut, special_lj, tstat_only, inum, nall, 300,
                       maxspecial, cell_size, gpu_split, screen);

  DPDMF.device->world_barrier();
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
      init_ok=DPDMF.init(ntypes, cutsq, host_a0, host_gamma, host_sigma,
                         host_cut, special_lj, tstat_only, inum, nall, 300,
                         maxspecial, cell_size, gpu_split, screen);

    DPDMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    DPDMF.estimate_gpu_overhead();
  return init_ok;
}

void dpd_gpu_clear() {
  DPDMF.clear();
}

int ** dpd_gpu_compute_n(const int ago, const int inum_full, const int nall,
                         double **host_x, int *host_type, double *sublo,
                         double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum, const double cpu_time, bool &success,
                         double **host_v, const double dtinvsqrt,
                         const int seed, const int timestep,
                         double *boxlo, double *prd) {
  return DPDMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success,
                       host_v, dtinvsqrt, seed, timestep, boxlo, prd);
}

void dpd_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success, tagint *tag,
                     double **host_v, const double dtinvsqrt,
                     const int seed, const int timestep,
                     const int nlocal, double *boxlo, double *prd) {
  DPDMF.compute(ago, inum_full, nall, host_x, host_type, ilist, numj,
                firstneigh, eflag, vflag, eatom, vatom, host_start, cpu_time, success,
                tag, host_v, dtinvsqrt, seed, timestep, nlocal, boxlo, prd);
}

void dpd_gpu_update_coeff(int ntypes, double **host_a0, double **host_gamma,
                          double **host_sigma, double **host_cut)
{
   DPDMF.update_coeff(ntypes,host_a0,host_gamma,host_sigma,host_cut);
}

double dpd_gpu_bytes() {
  return DPDMF.host_memory_usage();
}


