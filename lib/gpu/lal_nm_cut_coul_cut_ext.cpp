/***************************************************************************
                                 nm_cut_coul_cut_ext.cpp
                             -------------------

  Functions for LAMMPS access to nm/cut/coul/cut acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_nm_cut_coul_cut.h"
#include "lammps_gpu.h"

using namespace std;
using namespace LAMMPS_AL;

static NMCutCoulCut<PRECISION,ACC_PRECISION> NMCCMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------

namespace LAMMPS_GPU {
int nm_cut_coul_cut_gpu_init(const int ntypes, double **host_e0nm, double **host_nm,
                             double **host_r0n, double **host_r0m, double **host_nn,
                             double **host_mm, double **host_cut_ljsq, double **host_offset,
                             double **host_cutsq, double **host_cut_coulsq, double *special_lj,
                             double *special_coul, const int inum, const int nall,
                             const int max_nbors, const int maxspecial, const double cell_size,
                             int &gpu_mode, FILE *screen, const double host_qqrd2e)
{
  NMCCMF.clear();
  gpu_mode=NMCCMF.device->gpu_mode();
  int first_gpu=NMCCMF.device->first_device();
  int last_gpu=NMCCMF.device->last_device();
  int world_me=NMCCMF.device->world_me();
  int gpu_rank=NMCCMF.device->gpu_rank();
  int procs_per_gpu=NMCCMF.device->procs_per_gpu();

  NMCCMF.device->init_message(screen,"nm/cut/coul/cut",first_gpu,last_gpu);

  bool message=false;
  if (NMCCMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=NMCCMF.init(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm,
                        host_cut_ljsq, host_offset, host_cutsq, host_cut_coulsq, special_lj,
                        special_coul, inum, nall, max_nbors, maxspecial, cell_size, screen,
                        host_qqrd2e);

  NMCCMF.device->world_barrier();
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
      init_ok=NMCCMF.init(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm,
                          host_cut_ljsq, host_offset, host_cutsq, host_cut_coulsq, special_lj,
                          special_coul, inum, nall, max_nbors, maxspecial, cell_size, screen,
                          host_qqrd2e);

    NMCCMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    NMCCMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void nm_cut_coul_cut_gpu_reinit(const int ntypes, double **host_e0nm, double **host_nm,
                                double **host_r0n, double **host_r0m, double **host_nn,
                                double **host_mm, double **host_cut_ljsq, double **host_offset,
                                double **host_cutsq, double **host_cut_coulsq) {
  int world_me=NMCCMF.device->world_me();
  int gpu_rank=NMCCMF.device->gpu_rank();
  int procs_per_gpu=NMCCMF.device->procs_per_gpu();

  if (world_me==0)
    NMCCMF.reinit(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm, host_cut_ljsq,
                  host_offset, host_cutsq, host_cut_coulsq);
  NMCCMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      NMCCMF.reinit(ntypes, host_e0nm, host_nm, host_r0n, host_r0m, host_nn, host_mm, host_cut_ljsq,
                    host_offset, host_cutsq, host_cut_coulsq);
    NMCCMF.device->serialize_init();
  }
}

void nm_cut_coul_cut_gpu_clear() {
  NMCCMF.clear();
}

int **nm_cut_coul_cut_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                    double **host_x, int *host_type, double *sublo, double *subhi,
                                    tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                    const bool vflag, const bool eatom, const bool vatom,
                                    int **ilist, int **jnum, bool &success, double *host_q,
                                    double *boxlo, double *prd, int *periodicity)
{
  return NMCCMF.compute(ago, inum_full, nall, host_x, host_type, sublo, subhi, tag, nspecial,
                        special, eflag, vflag, eatom, vatom, ilist, jnum, success, host_q, boxlo,
                        prd, periodicity);
}

void nm_cut_coul_cut_gpu_compute(const int ago, const int inum_full, const int nall,
                                 double **host_x, int *host_type, int *ilist, int *numj,
                                 int **firstneigh, const bool eflag, const bool vflag,
                                 const bool eatom, const bool vatom, bool &success, double *host_q,
                                 const int nlocal, double *boxlo, double *prd)
{
  NMCCMF.compute(ago, inum_full, nall, host_x, host_type, ilist, numj, firstneigh, eflag, vflag,
                 eatom, vatom, success, host_q, nlocal, boxlo, prd);
}

double nm_cut_coul_cut_gpu_bytes() {
  return NMCCMF.host_memory_usage();
}

} // namespace LAMMPS_GPU
