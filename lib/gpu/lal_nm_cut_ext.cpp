/***************************************************************************
                                 nm_cut_ext.cpp
                             -------------------

  Functions for LAMMPS access to nm/cut acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_nm_cut.h"
#include "lammps_gpu.h"

using namespace std;
using namespace LAMMPS_AL;

static NMCut<PRECISION,ACC_PRECISION> NMCMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------

namespace LAMMPS_GPU {
int nm_cut_gpu_init(const int ntypes, double **host_e0nm, double **host_nm, double **host_r0n,
                    double **host_r0m, double **host_nn, double **host_mm, double **host_cutsq,
                    double **host_offset, double *special_lj, const int inum, const int nall,
                    const int max_nbors, const int maxspecial, const double cell_size,
                    int &gpu_mode, FILE *screen)
{
  NMCMF.clear();
  gpu_mode=NMCMF.device->gpu_mode();
  int first_gpu=NMCMF.device->first_device();
  int last_gpu=NMCMF.device->last_device();
  int world_me=NMCMF.device->world_me();
  int gpu_rank=NMCMF.device->gpu_rank();
  int procs_per_gpu=NMCMF.device->procs_per_gpu();

  NMCMF.device->init_message(screen,"nm/cut",first_gpu,last_gpu);

  bool message=false;
  if (NMCMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=NMCMF.init(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm, host_cutsq,
                       host_offset, special_lj, inum, nall, max_nbors, maxspecial, cell_size,
                       screen);

  NMCMF.device->world_barrier();
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
      init_ok=NMCMF.init(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm,
                         host_cutsq, host_offset, special_lj, inum, nall, max_nbors, maxspecial,
                         cell_size, screen);

    NMCMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    NMCMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void nm_cut_gpu_reinit(const int ntypes, double **host_e0nm, double **host_nm, double **host_r0n,
                       double **host_r0m, double **host_nn, double **host_mm, double **host_cutsq,
                       double **host_offset) {
  int world_me=NMCMF.device->world_me();
  int gpu_rank=NMCMF.device->gpu_rank();
  int procs_per_gpu=NMCMF.device->procs_per_gpu();

  if (world_me==0)
    NMCMF.reinit(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm, host_cutsq,
                 host_offset);
  NMCMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      NMCMF.reinit(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm, host_cutsq,
                   host_offset);
    NMCMF.device->serialize_init();
  }
}

void nm_cut_gpu_clear() {
  NMCMF.clear();
}

int **nm_cut_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int **ilist, int **jnum,
                           bool &success, double *prd, int *periodicity)
{
  return NMCMF.compute(ago, inum_full, nall, host_x, host_type, sublo, subhi, tag, nspecial,
                       special, eflag, vflag, eatom, vatom, ilist, jnum, success, prd,
                       periodicity);
}

void nm_cut_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                        int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                        const bool vflag, const bool eatom, const bool vatom, bool &success)
{
  NMCMF.compute(ago, inum_full, nall, host_x, host_type, ilist, numj, firstneigh, eflag, vflag,
                eatom, vatom, success);
}

double nm_cut_gpu_bytes() {
  return NMCMF.host_memory_usage();
}

} // namespace LAMMPS_GPU
