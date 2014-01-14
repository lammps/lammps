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
#include <math.h>

#include "lal_sw.h"

using namespace std;
using namespace LAMMPS_AL;

static SW<PRECISION,ACC_PRECISION> SWMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int sw_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors, 
                const double cell_size, int &gpu_mode, FILE *screen,
                int* host_map, const int nelements, int*** host_elem2param, const int nparams,
                const double* sw_epsilon, const double* sw_sigma,
                const double* sw_lambda, const double* sw_gamma,
                const double* sw_costheta, const double* sw_biga,
                const double* sw_bigb, const double* sw_powerp,
                const double* sw_powerq, const double* sw_cut, 
                const double* sw_cutsq) {
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
    init_ok=SWMF.init(ntypes, inum, nall, 300, cell_size, gpu_split, screen,
                      host_map, nelements, host_elem2param, nparams,
                      sw_epsilon, sw_sigma, sw_lambda, sw_gamma, sw_costheta, 
                      sw_biga, sw_bigb, sw_powerp, sw_powerq, sw_cut, sw_cutsq);

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
      init_ok=SWMF.init(ntypes, inum, nall, 300, cell_size, gpu_split, screen,
                        host_map, nelements, host_elem2param, nparams,
                        sw_epsilon, sw_sigma, sw_lambda, sw_gamma, sw_costheta, 
                        sw_biga, sw_bigb, sw_powerp, sw_powerq, sw_cut, 
                        sw_cutsq);

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
                        double *sublo, double *subhi, int *tag, int **nspecial,
                        int **special, const bool eflag, const bool vflag,
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


