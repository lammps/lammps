/***************************************************************************
                                 nm_cut_split_ext.cpp
                             -------------------

  Functions for LAMMPS access to nm/cut/split acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_nm_cut_split.h"
#include "lammps_gpu.h"

using namespace std;
using namespace LAMMPS_AL;

static NMCutSplit<PRECISION,ACC_PRECISION> NMCSMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------

namespace LAMMPS_GPU {
int nm_cut_split_gpu_init(const int ntypes, double **host_e0nm, double **host_nm, double **host_r0n,
                          double **host_r0m, double **host_nn, double **host_mm,
                          double **host_cutsq, double **host_offset, double **host_e0,
                          double **host_r0, double *special_lj, const int inum, const int nall,
                          const int max_nbors, const int maxspecial, const double cell_size,
                          int &gpu_mode, FILE *screen)
{
  NMCSMF.clear();
  gpu_mode=NMCSMF.device->gpu_mode();
  int first_gpu=NMCSMF.device->first_device();
  int last_gpu=NMCSMF.device->last_device();
  int world_me=NMCSMF.device->world_me();
  int gpu_rank=NMCSMF.device->gpu_rank();
  int procs_per_gpu=NMCSMF.device->procs_per_gpu();

  NMCSMF.device->init_message(screen,"nm/cut/split",first_gpu,last_gpu);

  bool message=false;
  if (NMCSMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=NMCSMF.init(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm,
                        host_cutsq, host_offset, host_e0, host_r0, special_lj, inum, nall,
                        max_nbors, maxspecial, cell_size, screen);

  NMCSMF.device->world_barrier();
  if (message)
    fprintf(screen,"Done.\n");

  for (int i=0; i<procs_per_gpu; i++) {
    if (message) {
      if (last_gpu-first_gpu==0)
        fprintf(screen,"Initializing Device %d on core %d...",first_gpu,i);
      else
        fprintf(screen,"Initializing Devices %d-%d on core %d...",first_gpu, last_gpu,i);
      fflush(screen);
    }
    if (gpu_rank==i && world_me!=0)
      init_ok=NMCSMF.init(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm,
                          host_cutsq, host_offset, host_e0, host_r0, special_lj, inum, nall,
                          max_nbors, maxspecial, cell_size, screen);

    NMCSMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    NMCSMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void nm_cut_split_gpu_reinit(const int ntypes, double **host_e0nm, double **host_nm,
                             double **host_r0n, double **host_r0m, double **host_nn,
                             double **host_mm, double **host_cutsq, double **host_offset,
                             double **host_e0, double **host_r0) {
  int world_me=NMCSMF.device->world_me();
  int gpu_rank=NMCSMF.device->gpu_rank();
  int procs_per_gpu=NMCSMF.device->procs_per_gpu();

  if (world_me==0)
    NMCSMF.reinit(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm, host_cutsq,
                  host_offset, host_e0, host_r0);
  NMCSMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      NMCSMF.reinit(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm, host_cutsq,
                    host_offset, host_e0, host_r0);
    NMCSMF.device->serialize_init();
  }
}

void nm_cut_split_gpu_clear() {
  NMCSMF.clear();
}

int **nm_cut_split_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                 double **host_x, int *host_type, double *sublo, double *subhi,
                                 tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                 int **jnum, bool &success, double *prd, int *periodicity)
{
  return NMCSMF.compute(ago, inum_full, nall, host_x, host_type, sublo, subhi, tag, nspecial,
                        special, eflag, vflag, eatom, vatom, ilist, jnum, success, prd,
                        periodicity);
}

void nm_cut_split_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, bool &success)
{
  NMCSMF.compute(ago, inum_full, nall, host_x, host_type, ilist, numj, firstneigh, eflag, vflag,
                 eatom, vatom, success);
}

double nm_cut_split_gpu_bytes() {
  return NMCSMF.host_memory_usage();
}

} // namespace LAMMPS_GPU
