/***************************************************************************
                                 born_ext.cpp
                             -------------------
                           Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to born acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_born.h"

using namespace std;
using namespace LAMMPS_AL;

static Born<PRECISION,ACC_PRECISION> BORNMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int born_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                  double **host_born1, double **host_born2, 
                  double **host_born3, double **host_a, double **host_c, 
                  double **host_d, double **sigma,      
                  double **offset, double *special_lj, const int inum,
                  const int nall, const int max_nbors,  const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen) {
  BORNMF.clear();
  gpu_mode=BORNMF.device->gpu_mode();
  double gpu_split=BORNMF.device->particle_split();
  int first_gpu=BORNMF.device->first_device();
  int last_gpu=BORNMF.device->last_device();
  int world_me=BORNMF.device->world_me();
  int gpu_rank=BORNMF.device->gpu_rank();
  int procs_per_gpu=BORNMF.device->procs_per_gpu();

  BORNMF.device->init_message(screen,"born",first_gpu,last_gpu);

  bool message=false;
  if (BORNMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=BORNMF.init(ntypes, cutsq, host_rhoinv, host_born1, host_born2, 
                        host_born3, host_a, host_c, host_d, sigma,
                        offset, special_lj, inum, nall, 300,
                        maxspecial, cell_size, gpu_split, screen);

  BORNMF.device->world_barrier();
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
      init_ok=BORNMF.init(ntypes, cutsq, host_rhoinv, host_born1, host_born2, 
                          host_born3, host_a, host_c, host_d, sigma, 
                          offset, special_lj, inum, nall, 300,
                          maxspecial, cell_size, gpu_split, screen);

    BORNMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    BORNMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void born_gpu_reinit(const int ntypes, double **host_rhoinv,
                     double **host_born1, double **host_born2,
                     double **host_born3, double **host_a, double **host_c,
                     double **host_d, double **offset) {
  int world_me=BORNMF.device->world_me();
  int gpu_rank=BORNMF.device->gpu_rank();
  int procs_per_gpu=BORNMF.device->procs_per_gpu();
  
  if (world_me==0)
    BORNMF.reinit(ntypes, host_rhoinv, host_born1, host_born2,
                  host_born3, host_a, host_c, host_d, offset);
  
  BORNMF.device->world_barrier();
  
  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      BORNMF.reinit(ntypes, host_rhoinv, host_born1, host_born2,
                    host_born3, host_a, host_c, host_d, offset);
    
    BORNMF.device->gpu_barrier();
  }
}

void born_gpu_clear() {
  BORNMF.clear(); 
}

int ** born_gpu_compute_n(const int ago, const int inum_full,
                          const int nall, double **host_x, int *host_type,
                          double *sublo, double *subhi, tagint *tag, int **nspecial,
                          tagint **special, const bool eflag, const bool vflag,
                          const bool eatom, const bool vatom, int &host_start,
                          int **ilist, int **jnum, const double cpu_time,
                          bool &success) {
  return BORNMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, ilist, jnum, cpu_time, success);
}  
			
void born_gpu_compute(const int ago, const int inum_full, const int nall,
                      double **host_x, int *host_type, int *ilist, int *numj,
                      int **firstneigh, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success) {
  BORNMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                 firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double born_gpu_bytes() {
  return BORNMF.host_memory_usage();
}


