/***************************************************************************
                               lal_ewald_ext.cpp
                             -------------------
                            W. Michael Brown (ORNL)
                            Axel Kohlmeyer (Temple)

  Functions for LAMMPS access to Ewald (reciprocal-space) acceleration routines

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : developers@lammps.org
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_ewald.h"

using namespace std;
using namespace LAMMPS_AL;

static EwaldGPU<PRECISION,ACC_PRECISION> EWALDMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
void ewald_gpu_init(const int nlocal, const int nall, FILE *screen,
                      int &success) {
  EWALDMF.clear(0.0);
  int first_gpu=EWALDMF.device->first_device();
  int last_gpu=EWALDMF.device->last_device();
  int world_me=EWALDMF.device->world_me();
  int gpu_rank=EWALDMF.device->gpu_rank();
  int procs_per_gpu=EWALDMF.device->procs_per_gpu();

  EWALDMF.device->init_message(screen,"ewald",first_gpu,last_gpu);

  bool message=false;
  if (EWALDMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  success=0;
  if (world_me==0)
    EWALDMF.init(nlocal,nall,screen,success);

  EWALDMF.device->world_barrier();
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
      EWALDMF.init(nlocal,nall,screen,success);

    EWALDMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");
}

// ---------------------------------------------------------------------------
// Upload the (constant per box) k-vectors and grid parameters
// ---------------------------------------------------------------------------
void ewald_gpu_setup(const int kmax, const int kcount, int *kxvecs,
                       int *kyvecs, int *kzvecs, double *ug, double **eg,
                       double **vg, double *unitk, int &success) {
  bool succ=true;
  EWALDMF.setup(kmax,kcount,kxvecs,kyvecs,kzvecs,ug,eg,vg,unitk,succ);
  success = succ ? 0 : -3;
}

// ---------------------------------------------------------------------------
// Compute the local (per-rank) structure factors on the device
// ---------------------------------------------------------------------------
int ewald_gpu_structure(const int ago, const int nlocal, const int nall,
                          double **host_x, int *host_type, double *host_q,
                          double *host_sfacrl, double *host_sfacim,
                          bool &success) {
  return EWALDMF.structure(ago,nlocal,nall,host_x,host_type,host_q,
                           host_sfacrl,host_sfacim,success);
}

// ---------------------------------------------------------------------------
// K-space field/force from the global structure factors
// ---------------------------------------------------------------------------
void ewald_gpu_compute(double *host_sfacrl_all, double *host_sfacim_all,
                         const double qscale, const int slabflag,
                         const int eflag_atom, const int vflag_atom,
                         double *host_eatom, double **host_vatom,
                         bool &success) {
  EWALDMF.compute_forces(host_sfacrl_all,host_sfacim_all,qscale,slabflag,
                         eflag_atom,vflag_atom,host_eatom,host_vatom,success);
}

void ewald_gpu_clear(const double cpu_time) {
  EWALDMF.clear(cpu_time);
}

double ewald_gpu_bytes() {
  return EWALDMF.host_memory_usage();
}
