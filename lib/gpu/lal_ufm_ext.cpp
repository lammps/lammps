/***************************************************************************
                                   ufm_ext.cpp
                           ------------------------------
                            Rodolfo Paula Leite (Unicamp/Brazil)
                            Maurice de Koning (Unicamp/Brazil)

  Functions for LAMMPS access to ufm acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : pl.rodolfo@gmail.com
                           dekoning@ifi.unicamp.br
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_ufm.h"

using namespace std;
using namespace LAMMPS_AL;

static UFM<PRECISION,ACC_PRECISION> UFMLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int ufml_gpu_init(const int ntypes, double **cutsq, double **host_uf1,
                 double **host_uf2, double **host_uf3, double **offset,
                 double *special_lj, const int inum, const int nall,
                 const int max_nbors,  const int maxspecial, const double cell_size,
                 int &gpu_mode, FILE *screen) {
  UFMLMF.clear();
  gpu_mode=UFMLMF.device->gpu_mode();
  double gpu_split=UFMLMF.device->particle_split();
  int first_gpu=UFMLMF.device->first_device();
  int last_gpu=UFMLMF.device->last_device();
  int world_me=UFMLMF.device->world_me();
  int gpu_rank=UFMLMF.device->gpu_rank();
  int procs_per_gpu=UFMLMF.device->procs_per_gpu();

  UFMLMF.device->init_message(screen,"ufm",first_gpu,last_gpu);

  bool message=false;
  if (UFMLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=UFMLMF.init(ntypes, cutsq, host_uf1, host_uf2, host_uf3,
                        offset, special_lj, inum, nall, max_nbors,
                        maxspecial, cell_size, gpu_split, screen);

  UFMLMF.device->world_barrier();
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
      init_ok=UFMLMF.init(ntypes, cutsq, host_uf1, host_uf2, host_uf3,
                         offset, special_lj, inum, nall, max_nbors, maxspecial,
                         cell_size, gpu_split, screen);

    UFMLMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    UFMLMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Copy updated coeffs from host to device
// ---------------------------------------------------------------------------
void ufml_gpu_reinit(const int ntypes, double **cutsq, double **host_uf1,
                    double **host_uf2, double **host_uf3, double **offset) {
  int world_me=UFMLMF.device->world_me();
  int gpu_rank=UFMLMF.device->gpu_rank();
  int procs_per_gpu=UFMLMF.device->procs_per_gpu();

  if (world_me==0)
    UFMLMF.reinit(ntypes, cutsq, host_uf1, host_uf2, host_uf3, offset);
  UFMLMF.device->world_barrier();

  for (int i=0; i<procs_per_gpu; i++) {
    if (gpu_rank==i && world_me!=0)
      UFMLMF.reinit(ntypes, cutsq, host_uf1, host_uf2, host_uf3, offset);
    UFMLMF.device->serialize_init();
  }
}

void ufml_gpu_clear() {
  UFMLMF.clear();
}

int ** ufml_gpu_compute_n(const int ago, const int inum_full,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum, const double cpu_time,
                        bool &success) {
  return UFMLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success);
}

void ufml_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success) {
  UFMLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
                firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double ufml_gpu_bytes() {
  return UFMLMF.host_memory_usage();
}


