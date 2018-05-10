/***************************************************************************
                                 lal_table.cpp
                             -------------------
                           Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to table acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_table.h"

using namespace std;
using namespace LAMMPS_AL;

static Table<PRECISION,ACC_PRECISION> TBMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int table_gpu_init(const int ntypes, double **cutsq, double ***table_coeffs,
                 double **table_data, double *special_lj, const int inum,
                 const int nall, const int max_nbors,  const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen,
                 int tabstyle, int ntables, int tablength) {
  TBMF.clear();
  gpu_mode=TBMF.device->gpu_mode();
  double gpu_split=TBMF.device->particle_split();
  int first_gpu=TBMF.device->first_device();
  int last_gpu=TBMF.device->last_device();
  int world_me=TBMF.device->world_me();
  int gpu_rank=TBMF.device->gpu_rank();
  int procs_per_gpu=TBMF.device->procs_per_gpu();

  TBMF.device->init_message(screen,"table",first_gpu,last_gpu);

  bool message=false;
  if (TBMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=TBMF.init(ntypes, cutsq, table_coeffs, table_data,
                      special_lj, inum, nall, 300, maxspecial, cell_size,
                      gpu_split, screen, tabstyle, ntables, tablength);

  TBMF.device->world_barrier();
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
      init_ok=TBMF.init(ntypes, cutsq, table_coeffs, table_data,
                      special_lj, inum, nall, 300, maxspecial, cell_size,
                      gpu_split, screen, tabstyle, ntables, tablength);

    TBMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    TBMF.estimate_gpu_overhead();
  return init_ok;
}

void table_gpu_clear() {
  TBMF.clear();
}

int ** table_gpu_compute_n(const int ago, const int inum_full,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum, const double cpu_time,
                        bool &success) {
  return TBMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success);
}

void table_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success) {
  TBMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double table_gpu_bytes() {
  return TBMF.host_memory_usage();
}


