/***************************************************************************
                                 lj96_ext.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Functions for LAMMPS access to lj96/cut acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_lj96.h"

using namespace std;
using namespace LAMMPS_AL;

static LJ96<PRECISION,ACC_PRECISION> LJ96MF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int lj96_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                  double **host_lj2, double **host_lj3, double **host_lj4,
                  double **offset, double *special_lj, const int inum,
                  const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen) {
  LJ96MF.clear();
  gpu_mode=LJ96MF.device->gpu_mode();
  double gpu_split=LJ96MF.device->particle_split();
  int first_gpu=LJ96MF.device->first_device();
  int last_gpu=LJ96MF.device->last_device();
  int world_me=LJ96MF.device->world_me();
  int gpu_rank=LJ96MF.device->gpu_rank();
  int procs_per_gpu=LJ96MF.device->procs_per_gpu();

  LJ96MF.device->init_message(screen,"lj96/cut",first_gpu,last_gpu);

  bool message=false;
  if (LJ96MF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=LJ96MF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3,
                        host_lj4, offset, special_lj, inum, nall, 300,
                        maxspecial, cell_size, gpu_split, screen);

  LJ96MF.device->world_barrier();
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
      init_ok=LJ96MF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                          offset, special_lj, inum,  nall, 300, maxspecial,
                          cell_size, gpu_split, screen);

    LJ96MF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    LJ96MF.estimate_gpu_overhead();
  return init_ok;
}

void lj96_gpu_clear() {
  LJ96MF.clear();
}

int** lj96_gpu_compute_n(const int ago, const int inum_full,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum, const double cpu_time,
                         bool &success) {
  return LJ96MF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, ilist, jnum, cpu_time, success);
}  
			
void lj96_gpu_compute(const int ago, const int inum_full, const int nall,
                      double **host_x, int *host_type, int *ilist, int *numj,
                      int **firstneigh, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success) {
  LJ96MF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,firstneigh,
                 eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double lj96_gpu_bytes() {
  return LJ96MF.host_memory_usage();
}


