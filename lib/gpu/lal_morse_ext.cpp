/***************************************************************************
                                  morse.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Functions for LAMMPS access to morse acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_morse.h"

using namespace std;
using namespace LAMMPS_AL;

static Morse<PRECISION,ACC_PRECISION> MORMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int mor_gpu_init(const int ntypes, double **cutsq,
                 double **host_lj1, double **host_lj2, double **host_lj3, 
                 double **host_lj4, double **offset, double *special_lj,
                 const int inum, const int nall, const int max_nbors, 
                 const int maxspecial, const double cell_size, int &gpu_mode,
                 FILE *screen) {
  MORMF.clear();
  gpu_mode=MORMF.device->gpu_mode();
  double gpu_split=MORMF.device->particle_split();
  int first_gpu=MORMF.device->first_device();
  int last_gpu=MORMF.device->last_device();
  int world_me=MORMF.device->world_me();
  int gpu_rank=MORMF.device->gpu_rank();
  int procs_per_gpu=MORMF.device->procs_per_gpu();

  MORMF.device->init_message(screen,"morse",first_gpu,last_gpu);

  bool message=false;
  if (MORMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=MORMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, 
                       host_lj4, offset, special_lj, inum, nall, 300,
                       maxspecial, cell_size, gpu_split, screen);

  MORMF.device->world_barrier();
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
    if (gpu_rank==i && world_me!=0)
      init_ok=MORMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                         offset, special_lj, inum, nall, 300, maxspecial,
                         cell_size, gpu_split, screen);

    MORMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    MORMF.estimate_gpu_overhead();
  return init_ok;
}

void mor_gpu_clear() {
  MORMF.clear();
}

int** mor_gpu_compute_n(const int ago, const int inum_full,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, int *tag, int **nspecial,
                        int **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum, const double cpu_time,
                        bool &success) {
  return MORMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success);
}  
			
void mor_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success) {
  MORMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double mor_gpu_bytes() {
  return MORMF.host_memory_usage();
}


