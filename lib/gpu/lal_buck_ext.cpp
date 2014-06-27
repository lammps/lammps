/***************************************************************************
                                 buck_ext.cpp
                             -------------------
                           Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to buck acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_buck.h"

using namespace std;
using namespace LAMMPS_AL;

static Buck<PRECISION,ACC_PRECISION> BUCKMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int buck_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                 double **host_buck1, double **host_buck2, 
                 double **host_a, double **host_c,       
                 double **offset, double *special_lj, const int inum,
                 const int nall, const int max_nbors,  const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen) {
  BUCKMF.clear();
  gpu_mode=BUCKMF.device->gpu_mode();
  double gpu_split=BUCKMF.device->particle_split();
  int first_gpu=BUCKMF.device->first_device();
  int last_gpu=BUCKMF.device->last_device();
  int world_me=BUCKMF.device->world_me();
  int gpu_rank=BUCKMF.device->gpu_rank();
  int procs_per_gpu=BUCKMF.device->procs_per_gpu();

  BUCKMF.device->init_message(screen,"buck",first_gpu,last_gpu);

  bool message=false;
  if (BUCKMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=BUCKMF.init(ntypes, cutsq, host_rhoinv, host_buck1, host_buck2, 
                       host_a, host_c, offset, special_lj, inum, nall, 300,
                       maxspecial, cell_size, gpu_split, screen);

  BUCKMF.device->world_barrier();
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
      init_ok=BUCKMF.init(ntypes, cutsq, host_rhoinv, host_buck1, host_buck2, 
                       host_a, host_c, offset, special_lj, inum, nall, 300,
                       maxspecial, cell_size, gpu_split, screen);

    BUCKMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    BUCKMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void buck_gpu_reinit(const int ntypes, double **cutsq, double **host_rhoinv,
                  double **host_buck1, double **host_buck2,
                  double **host_a, double **host_c, double **offset) {
  int world_me=BUCKMF.device->world_me();
  int gpu_rank=BUCKMF.device->gpu_rank();
  int procs_per_gpu=BUCKMF.device->procs_per_gpu();
  
  if (world_me==0)
    BUCKMF.reinit(ntypes, cutsq, host_rhoinv, host_buck1, host_buck2,
                  host_a, host_c, offset);
  
  BUCKMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      BUCKMF.reinit(ntypes, cutsq, host_rhoinv, host_buck1, host_buck2,
                    host_a, host_c, offset);
    
    BUCKMF.device->gpu_barrier();
  }
}

void buck_gpu_clear() {
  BUCKMF.clear(); 
}

int ** buck_gpu_compute_n(const int ago, const int inum_full,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum, const double cpu_time,
                        bool &success) {
  return BUCKMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success);
}  
			
void buck_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success) {
  BUCKMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double buck_gpu_bytes() {
  return BUCKMF.host_memory_usage();
}


