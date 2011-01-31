/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
 
/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#include <iostream>
#include <cassert>
#include <math.h>

#include "pppm_gpu_memory.h"

using namespace std;

static PPPMGPUMemory<PRECISION,ACC_PRECISION> PPPMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
bool pppm_gpu_init(const int nlocal, const int nall, FILE *screen) {
  PPPMF.clear();
  int first_gpu=PPPMF.device->first_device();
  int last_gpu=PPPMF.device->last_device();
  int world_me=PPPMF.device->world_me();
  int gpu_rank=PPPMF.device->gpu_rank();
  int procs_per_gpu=PPPMF.device->procs_per_gpu();

  PPPMF.device->init_message(screen,"pppm",first_gpu,last_gpu);

  bool message=false;
  if (PPPMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  if (world_me==0) {
    bool init_ok=PPPMF.init(nlocal,nall,screen);
    if (!init_ok)
      return false;
  }

  PPPMF.device->world_barrier();
  if (message)
    fprintf(screen,"Done.\n");

  for (int i=0; i<procs_per_gpu; i++) {
    if (message) {
      if (last_gpu-first_gpu==0)
        fprintf(screen,"Initializing GPU %d on core %d...",first_gpu,i);
      else
        fprintf(screen,"Initializing GPUs %d-%d on core %d...",first_gpu,
                last_gpu,i);
      fflush(screen);
    }
    if (gpu_rank==i && world_me!=0) {
      bool init_ok=PPPMF.init(nlocal,nall,screen);
      if (!init_ok)
        return false;
    }
    PPPMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");
  return true;
}

void pppm_gpu_clear() {
  PPPMF.clear();
}

void pppm_gpu_compute(const int ago, const int nlocal, const int nall,
                      double **host_x, int *host_type, bool &success,
                      double *host_q) {
  PPPMF.compute(ago,nlocal,nall,host_x,host_type,success,host_q);
}

double pppm_gpu_bytes() {
  return PPPMF.host_memory_usage();
}
