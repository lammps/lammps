/***************************************************************************
                               dipole_long_lj_ext.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to dipole/cut acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_dipole_long_lj.h"

using namespace std;
using namespace LAMMPS_AL;

static DipoleLongLJ<PRECISION,ACC_PRECISION> DPLJMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int dplj_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                 double **host_lj2, double **host_lj3, double **host_lj4,
                 double **offset, double *special_lj, const int inum,
                 const int nall, const int max_nbors, const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen,
                 double **host_cut_ljsq, const double host_cut_coulsq,
                 double *host_special_coul, const double qqrd2e,
                 const double g_ewald) {
  DPLJMF.clear();
  gpu_mode=DPLJMF.device->gpu_mode();
  double gpu_split=DPLJMF.device->particle_split();
  int first_gpu=DPLJMF.device->first_device();
  int last_gpu=DPLJMF.device->last_device();
  int world_me=DPLJMF.device->world_me();
  int gpu_rank=DPLJMF.device->gpu_rank();
  int procs_per_gpu=DPLJMF.device->procs_per_gpu();

  DPLJMF.device->init_message(screen,"lj/cut/dipole/long",first_gpu,last_gpu);

  bool message=false;
  if (DPLJMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=DPLJMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3,
                       host_lj4, offset, special_lj, inum, nall, max_nbors,
                       maxspecial, cell_size, gpu_split, screen, host_cut_ljsq,
                       host_cut_coulsq, host_special_coul, qqrd2e, g_ewald);

  DPLJMF.device->world_barrier();
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
      init_ok=DPLJMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                         offset, special_lj, inum, nall, max_nbors, maxspecial,
                         cell_size, gpu_split, screen, host_cut_ljsq,
                         host_cut_coulsq, host_special_coul, qqrd2e, g_ewald);

    DPLJMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    DPLJMF.estimate_gpu_overhead();
  return init_ok;
}

void dplj_gpu_clear() {
  DPLJMF.clear();
}

int** dplj_gpu_compute_n(const int ago, const int inum_full,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum, const double cpu_time,
                        bool &success, double *host_q, double **host_mu,
                        double *boxlo, double *prd) {
  return DPLJMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success,
                       host_q, host_mu, boxlo, prd);
}

void dplj_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success, double *host_q,
                     double **host_mu, const int nlocal, double *boxlo, double *prd) {
  DPLJMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,firstneigh,eflag,
                vflag,eatom,vatom,host_start,cpu_time,success,host_q,host_mu,
                nlocal,boxlo,prd);
}

double dplj_gpu_bytes() {
  return DPLJMF.host_memory_usage();
}


