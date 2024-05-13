/***************************************************************************
                               lj_gromacs_ext.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to lj/gromacs acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_lj_gromacs.h"

using namespace std;
using namespace LAMMPS_AL;

static LJGROMACS<PRECISION,ACC_PRECISION> LJGRMMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int ljgrm_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                   double **host_lj2, double **host_lj3, double **host_lj4,
                   double *special_lj, const int inum,
                   const int nall, const int max_nbors, const int maxspecial,
                   const double cell_size, int &gpu_mode, FILE *screen,
                   double **host_ljsw1, double **host_ljsw2, double **host_ljsw3,
                   double **host_ljsw4, double **host_ljsw5,
                   double **cut_inner, double **cut_inner_sq) {
  LJGRMMF.clear();
  gpu_mode=LJGRMMF.device->gpu_mode();
  double gpu_split=LJGRMMF.device->particle_split();
  int first_gpu=LJGRMMF.device->first_device();
  int last_gpu=LJGRMMF.device->last_device();
  int world_me=LJGRMMF.device->world_me();
  int gpu_rank=LJGRMMF.device->gpu_rank();
  int procs_per_gpu=LJGRMMF.device->procs_per_gpu();

  LJGRMMF.device->init_message(screen,"lj/gromacs",first_gpu,last_gpu);

  bool message=false;
  if (LJGRMMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    LJGRMMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                 special_lj, inum, nall, max_nbors, maxspecial, cell_size,
                 gpu_split, screen, host_ljsw1, host_ljsw2, host_ljsw3,
                 host_ljsw4, host_ljsw5, cut_inner, cut_inner_sq);

  LJGRMMF.device->world_barrier();
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
      init_ok=LJGRMMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                           special_lj, inum, nall, max_nbors, maxspecial, cell_size,
                           gpu_split, screen, host_ljsw1, host_ljsw2, host_ljsw3,
                           host_ljsw4, host_ljsw5, cut_inner, cut_inner_sq);

    LJGRMMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    LJGRMMF.estimate_gpu_overhead();
  return init_ok;
}

void ljgrm_gpu_clear() {
  LJGRMMF.clear();
}

int ** ljgrm_gpu_compute_n(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success) {
  return LJGRMMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                         subhi, tag, nspecial, special, eflag, vflag, eatom,
                         vatom, host_start, ilist, jnum, cpu_time, success);
}

void ljgrm_gpu_compute(const int ago, const int inum_full, const int nall,
                       double **host_x, int *host_type, int *ilist, int *numj,
                       int **firstneigh, const bool eflag, const bool vflag,
                       const bool eatom, const bool vatom, int &host_start,
                       const double cpu_time, bool &success) {
  LJGRMMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                  firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}


double ljgrm_gpu_bytes() {
  return LJGRMMF.host_memory_usage();
}


