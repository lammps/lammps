/***************************************************************************
                                 sw_ext.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Functions for LAMMPS access to sw acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Tue March 26, 2013
    email                : brownw@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_sw.h"

using namespace std;
using namespace LAMMPS_AL;

static SW<PRECISION,ACC_PRECISION> SWMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int sw_gpu_init(const int ntypes, const int inum, const int nall,
                const int max_nbors, const double cell_size, int &gpu_mode,
                FILE *screen, double **ncutsq, double **ncut, double **sigma,
                double **powerp, double **powerq, double **sigma_gamma,
                double **c1, double **c2, double **c3, double **c4,
                double **c5, double **c6, double ***lambda_epsilon,
                double ***costheta, const int *map, int ***e2param) {
  SWMF.clear();
  gpu_mode=SWMF.device->gpu_mode();
  double gpu_split=SWMF.device->particle_split();
  int first_gpu=SWMF.device->first_device();
  int last_gpu=SWMF.device->last_device();
  int world_me=SWMF.device->world_me();
  int gpu_rank=SWMF.device->gpu_rank();
  int procs_per_gpu=SWMF.device->procs_per_gpu();

  // disable host/device split for now
  if (gpu_split != 1.0)
    return -8;

  SWMF.device->init_message(screen,"sw/gpu",first_gpu,last_gpu);

  bool message=false;
  if (SWMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=SWMF.init(ntypes, inum, nall, max_nbors, cell_size, gpu_split,
                      screen, ncutsq, ncut, sigma, powerp, powerq,
                      sigma_gamma, c1, c2, c3, c4, c5, c6, lambda_epsilon,
                      costheta, map, e2param);

  SWMF.device->world_barrier();
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
      init_ok=SWMF.init(ntypes, inum, nall, max_nbors, cell_size, gpu_split,
                        screen, ncutsq, ncut, sigma, powerp, powerq,
                        sigma_gamma, c1, c2, c3, c4, c5, c6, lambda_epsilon,
                        costheta, map, e2param);

    SWMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    SWMF.estimate_gpu_overhead();
  return init_ok;
}

void sw_gpu_clear() {
  SWMF.clear();
}

int ** sw_gpu_compute_n(const int ago, const int inum_full,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum, const double cpu_time,
                        bool &success) {
  return SWMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success);
}

void sw_gpu_compute(const int ago, const int nlocal, const int nall,
                    const int nlist, double **host_x, int *host_type,
                    int *ilist, int *numj, int **firstneigh, const bool eflag,
                    const bool vflag, const bool eatom, const bool vatom,
                    int &host_start, const double cpu_time, bool &success) {
  SWMF.compute(ago,nlocal,nall,nlist,host_x,host_type,ilist,numj,
               firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double sw_gpu_bytes() {
  return SWMF.host_memory_usage();
}
