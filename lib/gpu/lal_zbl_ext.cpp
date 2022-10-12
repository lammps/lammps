/***************************************************************************
                                 zbl_ext.cpp
                             -------------------
                               Trung Dac Nguyen

  Functions for LAMMPS access to zbl acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndactrung@gmail.com
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_zbl.h"

using namespace std;
using namespace LAMMPS_AL;

static ZBL<PRECISION,ACC_PRECISION> ZBLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int zbl_gpu_init(const int ntypes, double **cutsq, double **host_sw1,
                 double **host_sw2, double **host_sw3, double **host_sw4, double **host_sw5,
                 double **host_d1a, double **host_d2a, double **host_d3a, double **host_d4a,
                 double **host_zze, double cut_globalsq, double cut_innersq, double cut_inner,
                 const int inum, const int nall, const int max_nbors,
                 const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen) {
  ZBLMF.clear();
  gpu_mode=ZBLMF.device->gpu_mode();
  double gpu_split=ZBLMF.device->particle_split();
  int first_gpu=ZBLMF.device->first_device();
  int last_gpu=ZBLMF.device->last_device();
  int world_me=ZBLMF.device->world_me();
  int gpu_rank=ZBLMF.device->gpu_rank();
  int procs_per_gpu=ZBLMF.device->procs_per_gpu();

  ZBLMF.device->init_message(screen,"zbl",first_gpu,last_gpu);

  bool message=false;
  if (ZBLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=ZBLMF.init(ntypes, cutsq, host_sw1, host_sw2, host_sw3, host_sw4,
                       host_sw5, host_d1a, host_d2a, host_d3a, host_d4a, host_zze,
                       cut_globalsq, cut_innersq, cut_inner,
                       inum, nall, max_nbors, maxspecial, cell_size, gpu_split, screen);

  ZBLMF.device->world_barrier();
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
      init_ok=ZBLMF.init(ntypes, cutsq, host_sw1, host_sw2, host_sw3, host_sw4,
                         host_sw5, host_d1a, host_d2a, host_d3a, host_d4a, host_zze,
                         cut_globalsq, cut_innersq, cut_inner,
                         inum, nall, max_nbors, maxspecial, cell_size, gpu_split, screen);

    ZBLMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    ZBLMF.estimate_gpu_overhead();
  return init_ok;
}

void zbl_gpu_clear() {
  ZBLMF.clear();
}

int ** zbl_gpu_compute_n(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success) {
  return ZBLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success);
}

void zbl_gpu_compute(const int ago, const int inum_full, const int nall,
                       double **host_x, int *host_type, int *ilist, int *numj,
                       int **firstneigh, const bool eflag, const bool vflag,
                       const bool eatom, const bool vatom, int &host_start,
                       const double cpu_time, bool &success) {
  ZBLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
               firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double zbl_gpu_bytes() {
  return ZBLMF.host_memory_usage();
}


